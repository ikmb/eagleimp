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

#ifndef FPGA_H
#define FPGA_H

#include <string>
#include <unordered_map>
#include <chrono>
#include <memory>
#include <thread>
#include <vector>
#include <map>

#include "Device.h"
#include "DeviceCategory.h"

#include "Buffer.h"
#include "Spinlock.h"

#ifdef USE_AD_FPGA

extern "C" {
#include <admxrc3.h>
}

#endif

namespace hybridsys {

using namespace std;

class FPGA : Device
{


public:

    enum class Result {
        Success,
        Cancelled,
        Timeout
    };

    enum class DMAChannelType {
        Disabled,
        AxiStreamMaster,
        AxiStreamSlave,
        AxiMmap
    };

    enum Window {
		WIN_USER = 0,
		WIN_USER_PREFETCH = 1,
		WIN_CONTROL_REGISTER = 2,
		WIN_GENERIC_BRIDGE_REGISTER = 3,
		WIN_MAX = 4
	};

    enum DMAEngineAddr {
    	DMAENGINE0 = 0x40,
    	DMAENGINE1 = 0x80,
    	DMAENGINE2 = 0xc0,
    	DMAENGINE3 = 0x100
    };

    enum BridgeRegisterOffset {
    	OFF_ABORT_STATUS       = 0x0,
		OFF_CLEANUP_FIFO       = 0x4,
		OFF_IRQ_STATUS_IRQ_ACK = 0x8,
		OFF_IRQ_ENABLE         = 0xC
    };

    const unsigned int REG_STATUS = 0x1F;

#ifdef USE_AD_FPGA
    using DeviceHandle = ADMXRC3_HANDLE;
    using BufferHandle = ADMXRC3_BUFFER_HANDLE;

    static constexpr DeviceHandle DeviceHandleInvalidValue = ADB3_HANDLE_INVALID_VALUE;
#else
    using DeviceHandle = void*;
    using BufferHandle = void*;
    const DeviceHandle DeviceHandleInvalidValue = nullptr;
#endif

    /**
     * @brief FPGA Instantiates the native API for the device specified by the parameters
     * @param bus The logical assigned PCI Express bus
     * @param slot The logical assigned PCI Express slot (also called device)
     */
    FPGA(int bus, int slot);
    virtual ~FPGA();

    FPGA(const FPGA&) = delete;
    FPGA(FPGA&& f);

    // Device interface
    int getBus() override;
    int getSlot() override;
    const string& getSerialNumber() override;
    int getIndex() const override;

    /**
     * @brief Queries the device whether it does contain a valid bitstream or not.
     * Most methods will throw a runtime error if the FPGA does not contain a valid
     * bitstream. Use this method to check whether an FPGA is ready for use.
     * @return A boolean value indicating whether a valid bitstream is currently loaded.
     */
    bool isConfigured() const;
    void declareDMAChannel(unsigned channel, DMAChannelType type);

    static constexpr size_t confBlockLength = (256/8);
    const void* getConfigurationData();

    Result writeDMA(const FPGABuffer &buffer, unsigned channel, chrono::milliseconds timeout = chrono::milliseconds::max(), size_t size = 0);
    Result readDMA(FPGABuffer &buffer, unsigned channel, chrono::milliseconds timeout = chrono::milliseconds::max());
    void writeReg(uint64_t address, void *source, size_t length, enum Window win = WIN_USER);
    void readReg(uint64_t address, void *target, size_t length, enum Window win = WIN_USER);

    /**
     * @brief Tries to cancel a running transaction on the specified DMA channel.
     * @param channel The DMA channel to cancel a transaction on
     *
     * This method tries to cancel a running transaction on the specified DMA channel.
     * It returns as soon as the underlying driver has been instructed to cancel the
     * transactions, so it may return before the transactions have actually been torn down.
     * It is safe to call this method if there is no running transaction.
     */
    void cancel(unsigned channel);

    // cancels all FPGA DMA transactions and prevents new transactions
    void abort();

    DeviceHandle& createThreadHandle();
    DeviceHandle& getHandle();
    void destroyThreadHandle();
    BufferHandle lockBuffer(const void *buffer, size_t size);
    void unlockBuffer(BufferHandle bufferHandle);

    static int findIndex(int bus, int slot);

    // lock the FPGA for exclusive usage
    void lock();
    // WARNING: Do not use the FPGA after unlock() unless you have locked it again with lock()!
    void unlock();

private:
    int bus;
    int slot;
    string serial;

    bool configInitialized;
    array<char, confBlockLength> confData;

    bool isTraced;
    int adIndex;
    DeviceHandle masterHandle;

    Spinlock transactionHandleLock;
    unordered_map<int, DeviceHandle> transactionHandles;
    Spinlock deviceHandleLock;
    unordered_map<pthread_t, DeviceHandle> deviceHandles;
    vector<bool> cancellationRequested;
    bool aborted;

    unsigned maxDMAChannels;
    vector<DMAChannelType> channelTypes;

    string lockfile;
    bool locked;
    int lockfd;
    bool lockfdwr;

//    string getLockDirectory() const;
//    string getLockFileName() const;
    void reset(bool force = false);
    void waitForLock(int inotify_fd, int watch, const string &lockfile);

    const string lockfileprefix = "hybridsys-fpga-";
    const string lockdir = "/opt/alphadata/lock";
};

}

#endif // FPGA_H
