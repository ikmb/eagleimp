/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt and Jan Christian Kässens,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of EagleImp.
 *
 *    EagleImp is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    EagleImp is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with EagleImp. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iomanip>
#include <sstream>
#include <system_error>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <future>
#include <thread>
#include <chrono>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <pthread.h>
#include <pwd.h>
#include <sys/inotify.h>
}

#include "DeviceCategory.h"
#include "ThreadUtils.h"

#include "FPGA.h"

namespace hybridsys {

using namespace std;

static const char *module_path = "/sys/bus/pci/drivers/adb3";

/* static */ int FPGA::findIndex(int bus, int slot) {

    int index = -1;

    stringstream device_path;
    device_path << module_path << "/";
    device_path << "0000:" << hex << setfill('0') << setw(2) << bus << ":" << hex << setfill('0') << setw(2) << slot << ".0";
    device_path << "/adb3i";

    DIR *instance_dir = opendir(device_path.str().c_str());
    if(!instance_dir) {
        throw system_error(errno, system_category(), "Cannot open specified device. Is the adb3 driver running?");
    }

    struct dirent *entry = nullptr;
    while((entry = readdir(instance_dir))) {
        if(entry->d_type == DT_DIR && strncmp(entry->d_name, "adb3i", sizeof("adb3i"))) {
            if(strlen(entry->d_name) == 6) {
                index = (entry->d_name[5] - 0x30) & 0xFF;
                break;
            }
        }
    }
    closedir(instance_dir);

    if(index == -1)
        throw runtime_error("Something's wrong, I couldn't find a suitable ADB3 driver instance");
    return index;
}

void FPGA::lock() {

    lockfile = lockdir + "/" + lockfileprefix + serial;

    while (!locked) {
        lockfdwr = false;
        mode_t m = umask(0); // read current file creation mask
        int fd = open(lockfile.c_str(), O_RDWR | O_CREAT, 0666);
        umask(m); // reset file creation mask to value before creating the lockfile
        if(fd < 0) { // failed to get write access to lockfile
            //cerr << "WARNING: Could not get write access to lockfile " << lockfile << ", but I may still be able to acquire a lock..." << endl;
            fd = open(lockfile.c_str(), O_RDONLY); // try to get read access
            if (fd < 0) // something clearly went wrong (no write access to lockfile dir?)
                throw system_error(errno, system_category(), "Could not open lockfile " + lockfile);
        } else
            lockfdwr = true;
        flock(fd, LOCK_EX);

        struct stat st0;
        fstat(fd, &st0);
        if(st0.st_nlink == 0) { // the file was deleted by another process -> try again
            //cerr << "WARNING: Lockfile was deleted. Try to create a new one." << endl;
            close(fd);
            continue;
        }

        locked = true;
        lockfd = fd;
    }

    // Write some information to lockfile
    if (lockfdwr) {
        int ftr = ftruncate(lockfd, 0); // clear lockfile (don't care if this fails)

        stringstream buf;
        time_t epoch = time(nullptr);
        tm local = *localtime(&epoch);
        buf << "Since: " << put_time(&local, "%c") << endl;

        auto uid = geteuid();
        auto pw = getpwuid(uid);
        buf << "By: " << pw->pw_name << " (" << pw->pw_gecos << ")" << endl;

        buf << "Pid: " << getpid() << endl;
        char *cmdline = new char[4096]; // arbitrary
        int cmdline_fd = open("/proc/self/cmdline", O_RDONLY);
        if(cmdline_fd != -1) {
            cmdline[0] = 0;
            cmdline[4095] = 0;
            if(read(cmdline_fd, cmdline, 4096) == -1)
                throw system_error(errno, system_category(), "Could not read process details from /proc/self");
            buf << "Cmd: " << cmdline << endl;
        }
        delete[] cmdline;

        int ret = write(lockfd, buf.str().c_str(), buf.str().length()); // write to lockfile
        if (ftr || ret <= 0)
            cerr << "WARNING! Writing to lockfile failed with status " << ftr << "/" << ret << "." << endl;

        // DON'T DO THAT! THIS WILL RELEASE THE LOCK!
        // close(fd);
    }

#ifdef USE_AD_FPGA
    auto status = ADMXRC3_OpenEx(adIndex, false, 0, &masterHandle);
    switch(status) {
    case ADMXRC3_SUCCESS:
        break;
    case ADMXRC3_DEVICE_NOT_FOUND:
        // This device does not have a Tandem stage 2 bitstream loaded
        // That's okay, but we can't use the ADMXRC3 API, yet.
        break;
    default:
        throw system_error(status, FPGACategory());
    }
#endif

    reset(false);

}

void FPGA::unlock() {
    if (locked) {
        if (lockfdwr) {
            int ftr = ftruncate(lockfd, 0); // clear the contents of the lockfile
            if (ftr)
                cerr << "WARNING! Clearing FPGA lock file returned status " << ftr << "." << endl;
        }
        // important to remove the lockfile only while I have the exclusive lock on it!!!
        int unl = unlink(lockfile.c_str());
        if (unl)
            cerr << "WARNING! Failed to remove lock file." << endl;
        close(lockfd);
        lockfd = 0;
        locked = false;
        lockfdwr = false;
    }
}

#ifdef USE_AD_FPGA
static bool isDebuggerAttached() {
    ifstream f("/proc/self/status");
    string line;
    int tracerPid = 0;
    static const char token[] = "TracerPid:";

    while(getline(f, line)) {
        size_t pos = line.find(token);
        if(pos != string::npos) {
            tracerPid = atoi(line.data() + pos + sizeof(token));
        }
    }
    return tracerPid > 0;
}
#endif

FPGA::FPGA(int bus_, int slot_)
    :
      bus(bus_),
      slot(slot_),
      serial(""),
      configInitialized(false),
      confData(),
      adIndex(-1),
      masterHandle(FPGA::DeviceHandleInvalidValue),
      transactionHandleLock(),
      transactionHandles(),
      deviceHandleLock(),
      deviceHandles(),
      cancellationRequested(),
      aborted(false),
      maxDMAChannels(0),
      channelTypes(),
      locked(false),
      lockfd(0)
{
#ifdef USE_AD_FPGA
    ADMXRC3_STATUS status;

    adIndex = findIndex(bus_, slot_);

    status = ADMXRC3_OpenEx(adIndex, true, 0, &masterHandle);
    switch(status) {
    case ADMXRC3_SUCCESS:

        ADMXRC3_CARD_INFOEX info;
        memset(&info, 0, sizeof(info));
        ADMXRC3_GetCardInfoEx(masterHandle, &info);

        // not exactly precise but the only thing we know for sure
        if(info.NumDmaChannel > 0) {
            serial = to_string(info.SerialNumber);
            maxDMAChannels = info.NumDmaChannel;
            cancellationRequested.resize(maxDMAChannels);
            fill(begin(cancellationRequested), end(cancellationRequested), false);
            channelTypes.resize(maxDMAChannels);
            fill(begin(channelTypes), end(channelTypes), DMAChannelType::Disabled);
        } else {
            ADMXRC3_Close(masterHandle);
            masterHandle = ADB3_HANDLE_INVALID_VALUE;
        }
        break;
    case ADMXRC3_DEVICE_NOT_FOUND:
        // This device does not have a Tandem stage 2 bitstream loaded
        // That's okay, but we can't use the ADMXRC3 API, yet.
        break;
    default:
        throw system_error(status, FPGACategory());
    }

    /* Debuggers may emit signals that disturb some of the DMA calls.
     * They may return with CANCELLED statuses although no API cancellation
     * was requested. Knowing of this potential disturbance may reduce
     * false positives when the cause of cancellation is determined. */
    this->isTraced = isDebuggerAttached();
#endif
}

void FPGA::reset(bool force) {
    if(!locked)
        throw runtime_error("reset: This operation requires an FPGA lock");

    // check status if a reset is necessary
    createThreadHandle();
    unsigned char fpga_status; // one word
    readReg(REG_STATUS, &fpga_status, 1);

    if (fpga_status || force) { // some unclean status != 0 or forced

        if (!force)
            cerr << "WARNING: FPGA initialization detected unclean FPGA status (0x" << hex << setfill('0') << setw(2) << (unsigned int) fpga_status << dec << "). Trying a reset first." << endl;
        else
            cerr << "Forcing FPGA reset (curr status 0x" << hex << setfill('0') << setw(2) << (unsigned int) fpga_status << dec << ")." << endl;

        uint32_t init_reset = 1;
        writeReg(0, &init_reset, sizeof(init_reset), WIN_USER); // initialize reset
        init_reset = 0;
        writeReg(0, &init_reset, sizeof(init_reset), WIN_USER); // deassert reset flag

        // wait a bit
        this_thread::sleep_for(chrono::milliseconds(2000));

        // flush FIFOs
        uint32_t dummyword = 0;
        writeReg(DMAENGINE0+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);
        writeReg(DMAENGINE1+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);
        writeReg(DMAENGINE2+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);
        writeReg(DMAENGINE3+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);

        // wait again (but not so long)
        this_thread::sleep_for(chrono::milliseconds(500));

        readReg(REG_STATUS, &fpga_status, 1);
        if (fpga_status) {
            cerr << "WARNING: Unclean FPGA status after reset and FIFO flush (0x" << hex << setfill('0') << setw(2) << (unsigned int) fpga_status << dec << "). Continuing anyway." << endl;
        }

    } else {

        // flush FIFOs
        uint32_t dummyword = 0;
        writeReg(DMAENGINE0+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);
        writeReg(DMAENGINE1+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);
        writeReg(DMAENGINE2+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);
        writeReg(DMAENGINE3+OFF_CLEANUP_FIFO, &dummyword, sizeof(dummyword), WIN_GENERIC_BRIDGE_REGISTER);

    }
}

FPGA::FPGA(FPGA&& f)
    :
      bus(f.bus),
      slot(f.slot),
      serial(f.serial),
      configInitialized(f.configInitialized),
      confData(f.confData),
      isTraced(f.isTraced),
      adIndex(f.adIndex),
      masterHandle(f.masterHandle),
      transactionHandles(f.transactionHandles),
      deviceHandles(f.deviceHandles),
      cancellationRequested(f.cancellationRequested),
      aborted(f.aborted),
      maxDMAChannels(f.maxDMAChannels),
      channelTypes(f.channelTypes),
      locked(f.locked),
      lockfd(f.lockfd)
{
    f.masterHandle = FPGA::DeviceHandleInvalidValue;
    f.adIndex = -1;
}


FPGA::~FPGA() {
#ifdef USE_AD_FPGA
    if(masterHandle != FPGA::DeviceHandleInvalidValue) {
        bool transactionsActive = false;

        // initiate cancellations for active transactions, if any
        transactionHandleLock.lock();

        //cerr << "~FPGA(): There are " << transactionHandles.size() << " transaction handles." << endl;

        for(auto trans : transactionHandles) {
        	//cerr << " " << trans.second << flush;
            if(trans.second != ADB3_HANDLE_INVALID_VALUE) {
                transactionsActive = true;
                ADMXRC3_Cancel(trans.second);
            }
        }

        // Wait for cancelled transactions to actually be torn down.
        // Their respective threads will set the handle to an invalid value when done.
        // Probably more efficient using condition variables but this destructor is not time-critical.
        while(transactionsActive) {
            transactionsActive = false;
            for(auto trans : transactionHandles) {
                if(trans.second != ADB3_HANDLE_INVALID_VALUE) {
                    // still active, wait for cancellation
                    transactionsActive = true;
                }
            }
            transactionHandleLock.unlock();
            this_thread::yield();
            transactionHandleLock.lock();
        }
        transactionHandleLock.unlock();

        // No more transactions running, can now start closing any open file descriptors
        ADMXRC3_Close(masterHandle);

        for(auto entry : deviceHandles) {
            if(entry.second != ADB3_HANDLE_INVALID_VALUE)
                ADMXRC3_Close(entry.second);
        }
    }

    unlock();
#endif
}

int FPGA::getBus() { return bus; }

int FPGA::getSlot() { return slot; }

const string &FPGA::getSerialNumber() {
#ifdef USE_AD_FPGA
    if(masterHandle == FPGA::DeviceHandleInvalidValue)
        throw system_error(ADMXRC3_BAD_DRIVER, FPGACategory());
    return serial;
#else
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

int FPGA::getIndex() const {
#ifdef USE_AD_FPGA
    if(masterHandle == FPGA::DeviceHandleInvalidValue)
        throw system_error(ADMXRC3_BAD_DRIVER, FPGACategory());

    return adIndex;
#else
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

bool FPGA::isConfigured() const {
    return (masterHandle != FPGA::DeviceHandleInvalidValue);
}

FPGA::DeviceHandle& FPGA::createThreadHandle() {
#ifdef USE_AD_FPGA
    lock_guard<Spinlock> guard(deviceHandleLock);

    auto it = deviceHandles.find(pthread_self());
    if(it != end(deviceHandles))
        //throw system_error(ADMXRC3_INVALID_HANDLE, FPGACategory(), "This thread has already been initialized");
    	return it->second;

    FPGA::DeviceHandle newHandle (FPGA::DeviceHandleInvalidValue);
    ADMXRC3_STATUS status = ADMXRC3_OpenEx(adIndex, false, 0, &newHandle);

    if(status != ADMXRC3_SUCCESS)
        throw system_error(status, FPGACategory(), "Failed to open device");

    auto result = deviceHandles.emplace(make_pair(pthread_self(), newHandle));

    return (result.first)->second;
#else
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

FPGA::DeviceHandle& FPGA::getHandle() {
#ifdef USE_AD_FPGA
    lock_guard<Spinlock> guard(deviceHandleLock);

    unordered_map<pthread_t, FPGA::DeviceHandle>::iterator it;

    it = deviceHandles.find(pthread_self());
    if(it == end(deviceHandles)) {
        throw system_error(FPGA::DeviceHandleInvalidValue, FPGACategory(), "This thread must be initialized with the API before use");
    }

    return it->second;
#else
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

void FPGA::destroyThreadHandle() {
#ifdef USE_AD_FPGA
    lock_guard<Spinlock> guard(deviceHandleLock);

    auto it = deviceHandles.find(pthread_self());
    if(it != end(deviceHandles)) {
        ADMXRC3_Close(it->second);
        deviceHandles.erase(it);
    } else {
        throw system_error(FPGA::DeviceHandleInvalidValue, FPGACategory(), "This thread must be initialized with the API before use");
    }
#else
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

FPGA::BufferHandle FPGA::lockBuffer(const void *buffer, size_t size)
{
#ifdef USE_AD_FPGA
    BufferHandle bufferHandle;
    ADMXRC3_STATUS status = ADMXRC3_Lock(masterHandle, buffer, size, &bufferHandle);
    if(status == ADMXRC3_SUCCESS)
        return bufferHandle;
    throw system_error(status, FPGACategory(), "Failed locking buffer");
#else
    (void) buffer;
    (void) size;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

void FPGA::unlockBuffer(BufferHandle bufferHandle) {
#ifdef USE_AD_FPGA
    ADMXRC3_Unlock(masterHandle, bufferHandle);
#else
    (void) bufferHandle;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

FPGA::Result FPGA::writeDMA(const FPGABuffer &buffer, unsigned channel, chrono::milliseconds timeout, size_t size) {
#ifdef USE_AD_FPGA
    if(!locked)
        throw runtime_error("writeDMA: This operation requires an FPGA lock");

    if (aborted)
        return Result::Cancelled;

    Result ret = Result::Success;

    if(channel > channelTypes.size())
        throw system_error(ADMXRC3_INVALID_INDEX, FPGACategory(), "Undeclared DMA channel");

    if(channelTypes[channel] != DMAChannelType::AxiStreamMaster)
        throw system_error(ADMXRC3_INVALID_INDEX, FPGACategory(), "This DMA channel is not suitable for stream writing.");



    packaged_task<ADMXRC3_STATUS(void)> task([&]() {
        {
			stringstream ss;
			ss << "FPGA-" << adIndex << "-" << channel << "-DMA";
			ThreadUtils::setThreadName(ss.str());
        }

        ADMXRC3_STATUS finish_status = ADMXRC3_INTERNAL_ERROR;
        FPGA::DeviceHandle &handle = createThreadHandle();
        transactionHandleLock.lock();
        transactionHandles[channel] = handle;
        transactionHandleLock.unlock();

        try {
            if(size == 0)
                size = buffer.getSize();

            /*
            ADMXRC3_STATUS status = ADMXRC3_StartReadDMALockedEx(handle,
                                                                 NULL, // ticket, must be NULL for non-overlapping transfers
                                                                 channel,
                                                                 0, // flags
                                                                 buffer.getHandle(), // buffer handle
                                                                 0, // offset
                                                                 size,
                                                                 0);
                                                                 */
            ADMXRC3_STATUS status = ADMXRC3_StartWriteDMAEx(handle, NULL, channel, 0, buffer.getData(), size, 0);

            if(status != ADMXRC3_PENDING) {
                stringstream ss;
                ss << "ADAPI for FPGA " << adIndex << " failed while trying to write to channel " << channel << ": " << ADMXRC3_GetStatusString(status, false);
                throw system_error(status, FPGACategory(), ss.str());
            }

            finish_status = ADMXRC3_CANCELLED;
            do {
                finish_status = ADMXRC3_FinishDMA(handle, NULL, 1);

            } while(finish_status == ADMXRC3_CANCELLED && !cancellationRequested[channel] && isTraced);


            switch(finish_status) {
            case ADMXRC3_CANCELLED:
            case ADMXRC3_SUCCESS:
                break;
            default:
                throw system_error(finish_status, FPGACategory(), "Could not finish non-blocking DMA transaction");
            }

            destroyThreadHandle();

            transactionHandleLock.lock();
            transactionHandles[channel] = FPGA::DeviceHandleInvalidValue;
            transactionHandleLock.unlock();
        } catch(...) {
            destroyThreadHandle();

            transactionHandleLock.lock();
            transactionHandles[channel] = FPGA::DeviceHandleInvalidValue;
            transactionHandleLock.unlock();
            throw;
        }

        return finish_status;
    });

    future<ADMXRC3_STATUS> f = task.get_future();
    thread th (move(task));

    if(f.wait_for(timeout) == future_status::timeout) {
        cancel(channel);
        ret = Result::Timeout;
    }

    ADMXRC3_STATUS status = f.get();   // here, exceptions are thrown that happen within the packaged task
    th.join();

    if(ret != Result::Timeout && (status == ADMXRC3_CANCELLED || status == ADMXRC3_NONBLOCK_IDLE))
        ret = Result::Cancelled;

//    if(status != ADMXRC3_SUCCESS && ret != Result::Cancelled) {
//        stringstream ss;
//        ss << "ADAPI for FPGA " << adIndex << " failed while trying to write to channel " << channel << ": " << ADMXRC3_GetStatusString(status, false);
//
//        throw system_error(status, FPGACategory(), ss.str());
//    }

    cancellationRequested[channel] = false; // use atomic bool?

    return ret;
#else
    (void) buffer; (void) channel; (void) timeout; (void) size;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

void FPGA::declareDMAChannel(unsigned channel, DMAChannelType type) {
#ifdef USE_AD_FPGA
    if(channel > maxDMAChannels) {
        stringstream ss;
        ss << "The FPGA core does only support " << maxDMAChannels << " DMA channels. You cannot declare more than that.";
        throw system_error(ADMXRC3_INVALID_INDEX, FPGACategory(), ss.str());
    }

    channelTypes[channel] = type;
#else
    (void) channel; (void) type;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

const void* FPGA::getConfigurationData() {
    if(!configInitialized) {
    	confData.fill(0x00);
    	readReg(0, confData.data(), confData.size());
    	configInitialized = true;
    }
    return confData.data();
}

void FPGA::cancel(unsigned channel) {
#ifdef USE_AD_FPGA
    lock_guard<Spinlock> l(transactionHandleLock);

    auto it = transactionHandles.find(channel);
    if(it != end(transactionHandles)) {
        FPGA::DeviceHandle& handle = it->second;

        // In case there is no transaction present on this channel,
        // there will either be no entry in transactionHandles or a
        // ADMXRC3_HANDLE_INVALID_VALUE. It's okay to call Cancel with
        // that value, we're ignoring the return value anyway.
        cancellationRequested[channel] = true;
        ADMXRC3_Cancel(handle);
    }
#else
    (void) channel;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

void FPGA::abort(){
    aborted = true;
#ifdef USE_AD_FPGA
    lock_guard<Spinlock> l(transactionHandleLock);

    for (auto it = transactionHandles.begin(); it != transactionHandles.end(); ++it) {
        cancellationRequested[it->first] = true;
        ADMXRC3_Cancel(it->second);
    }
#else
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

FPGA::Result FPGA::readDMA(FPGABuffer &buffer, unsigned channel, chrono::milliseconds timeout) {
#ifdef USE_AD_FPGA
    if(!locked)
        throw runtime_error("readDMA: This operation requires an FPGA lock");

    if (aborted)
        return Result::Cancelled;

    Result ret = Result::Success;

    if(channel > channelTypes.size())
        throw system_error(ADMXRC3_INVALID_INDEX, FPGACategory(), "Undeclared DMA channel");

    if(channelTypes[channel] != DMAChannelType::AxiStreamSlave)
        throw system_error(ADMXRC3_INVALID_INDEX, FPGACategory(), "This DMA channel is not suitable for stream reading.");



    packaged_task<ADMXRC3_STATUS(void)> task([&]() {
        {
        	stringstream ss;
        	ss << "FPGA-" << adIndex << "-" << channel << "-DMA";
        	ThreadUtils::setThreadName(ss.str());
        }

        ADMXRC3_STATUS finish_status = ADMXRC3_INTERNAL_ERROR;
        FPGA::DeviceHandle &handle = createThreadHandle();
        transactionHandleLock.lock();
        transactionHandles[channel] = handle;
        transactionHandleLock.unlock();


        try {
            ADMXRC3_STATUS status = ADMXRC3_StartReadDMAEx(handle, NULL, channel, 0, buffer.getData(), buffer.getSize(), 0);

            if(status != ADMXRC3_PENDING) {
                stringstream ss;
                ss << "ADAPI for FPGA " << adIndex << " failed while trying to read from channel " << channel << ": " << ADMXRC3_GetStatusString(status, false);
                throw system_error(status, FPGACategory(), ss.str());
            }

            finish_status = ADMXRC3_CANCELLED;
            do {
                finish_status = ADMXRC3_FinishDMA(handle, NULL, 1);

                // tracing may send spurious cancellation signals
            } while(finish_status == ADMXRC3_CANCELLED && !cancellationRequested[channel] && isTraced);

            switch(finish_status) {
            case ADMXRC3_CANCELLED:
            case ADMXRC3_SUCCESS:
                break;
            default:
                throw system_error(finish_status, FPGACategory(), "Could not finish non-blocking DMA transaction");
            }

            destroyThreadHandle();
            transactionHandleLock.lock();
            transactionHandles[channel] = FPGA::DeviceHandleInvalidValue;
            transactionHandleLock.unlock();
        } catch(...) {
            destroyThreadHandle();
            transactionHandleLock.lock();
            transactionHandles[channel] = FPGA::DeviceHandleInvalidValue;
            transactionHandleLock.unlock();
            throw;
        }

        return finish_status;
    });

    future<ADMXRC3_STATUS> f = task.get_future();
    thread th (move(task));

    if(f.wait_for(timeout) == future_status::timeout) {
        cancel(channel);
        ret = Result::Timeout;
    }

    ADMXRC3_STATUS status = f.get();   // here, exceptions are thrown that happen within the packaged task
    th.join();

    if(ret != Result::Timeout && (status == ADMXRC3_CANCELLED || status == ADMXRC3_NONBLOCK_IDLE))
        ret = Result::Cancelled;

//    if(status != ADMXRC3_SUCCESS && ret != Result::Cancelled) {
//        stringstream ss;
//        ss << "ADAPI for FPGA " << adIndex << " failed while trying to read from channel " << channel << ": " << ADMXRC3_GetStatusString(status, false);
//
//        throw system_error(status, FPGACategory(), ss.str());
//    }

    cancellationRequested[channel] = false; // use atomic bool?

    return ret;
#else
    (void) buffer; (void) channel; (void) timeout;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

void FPGA::writeReg(uint64_t address, void *source, size_t length, enum Window win) {
#ifdef USE_AD_FPGA
    FPGA::DeviceHandle& handle = getHandle();
    ADMXRC3_STATUS status;
    status = ADMXRC3_Write(handle,
                 win,
                 0,
                 address,
                 length,
                 source);

    if (status != ADMXRC3_SUCCESS) {
        stringstream ss;
        ss << "ADAPI for FPGA " << this->adIndex << " failed while trying to write to register in window " << win << " at 0x" << hex << address << dec << ": " << ADMXRC3_GetStatusString(status, false);
        throw system_error(status, FPGACategory(), ss.str());
    }
#else
    (void) address; (void) source; (void) length; (void) win;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

void FPGA::readReg(uint64_t address, void *target, size_t length, enum Window win) {
#ifdef USE_AD_FPGA
    FPGA::DeviceHandle& handle = getHandle();
    ADMXRC3_STATUS status;
    status = ADMXRC3_Read(handle,
                 win,
                 0,
                 address,
                 length,
                 target);

    if (status != ADMXRC3_SUCCESS) {
        stringstream ss;
        ss << "ADAPI for FPGA " << this->adIndex << " failed while trying to read from register in window " << win << " at 0x" << hex << address << dec << ": " << ADMXRC3_GetStatusString(status, false);
        throw system_error(status, FPGACategory(), ss.str());
    }
#else
    (void) address; (void) target; (void) length; (void) win;
    throw runtime_error("Alpha Data FPGA API not available!");
#endif
}

} // namespace hybridsys
