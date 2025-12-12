//
// Created by marijn on 10/28/23.
//

#ifndef STD_LIB_STRING_H
#define STD_LIB_STRING_H

#include "defs.h"

#define strtok_next() strtok(NULL, NULL)


// string.asm
uint64_t strlen(const char* str);
char* strcpy(char* dst, const char* src);
char* strcat(char* dst, const char* src);
char* strtok(char* str, const char* delims);
// TODO: sse4: https://github.com/aklomp/sse-strings/blob/master/src/strcmpeq_sse4.asm


#endif //STD_LIB_STRING_H
