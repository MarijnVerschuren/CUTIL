//
// Created by marijn on 3/3/23.
//

#ifndef STD_LIB_SYS_H
#define STD_LIB_SYS_H

#include "defs.h"


// TODO: refine these
// https://wiki.osdev.org/CPUID
// https://en.wikipedia.org/wiki/CPUID
// https://c9x.me/x86/html/file_module_x86_id_45.html
enum CPUID_LEAFS {
	VENDOR_ID =										0x00000000,
	CPU_INFO =										0x00000001,
	TLB_INFO =										0x00000002,	// TODO
	SERIAL_NUMBER =									0x00000003,
	DETERMINISTIC_CACHE_INFO =						0x00000004,	// TODO
	MWAIT_MONITOR_INFO =							0x00000005,	// TODO
	THERMAL_INFO =									0x00000006,	// TODO
	EXTENDED_FEATURE_INFO_0 =						0x00000007,	// TODO
	// TODO: find purpose of results between 0x7 and 0x16
	END_STANDARD_RANGE =							0x00000016,	// the result of this is repeated when increasing eax

	EXTENDED_LEAF_INFO =							0x80000000,	// TODO
	EXTENDED_FEATURE_INFO_1 =						0x80000001,	// TODO
	BRAND_STRING0 =									0x80000002,
	BRAND_STRING1 =									0x80000003,
	BRAND_STRING2 =									0x80000004,
	L1_CACHE_INFO =									0x80000005,	// TODO
	L2_CACHE_INFO =									0x80000006,	// TODO
	POWER_MANAGEMENT_INFO =							0x80000007,	// TODO
	ADDRESS_SIZES =									0x80000008,	// TODO
	// TODO: find purpose of results between 0x8 and 0x1F
	ENCRYPTED_MEMORY_INFO =							0x8000001F,	// TODO
	EXTENDED_FEATURE_INFO_2 =						0x80000021,	// TODO

	AMD_EASTER_EGG =								0x8FFFFFFF,	// TODO
};

extern void cpuid(uint64_t leaf, uint32_t* eax, uint32_t* ebx, uint32_t* ecx, uint32_t* edx);


#endif //STD_LIB_SYS_H
