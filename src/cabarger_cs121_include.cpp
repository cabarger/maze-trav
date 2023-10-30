//!
//! cabarger_cs121_include.cpp
//!
//! CS 121.Bolden
//! Clang 16.0.6
//! Caleb Barger
//! 09/25/23 
//! Linux x86_64
//! barg8397@vandals.uidaho.edu
//!
//! Implementation file for cs121_include.h
//!

#include "cabarger_cs121_include.h"
#include "stdlib.h"
#include <string.h> // memcmp
#include <stdio.h> // FILE*
#include <time.h> 

void* arenaAllocInner(Arena* a, const u64 requested_space) {
  if (a->offset + requested_space > a->cap) 
    AssertMessage("out of space");
  void* result = a->base_ptr + a->offset;  
  a->offset += requested_space; 
  return result; 
}

Arena arenaInit(const u64 cap) {
  return Arena{
    .base_ptr = (u8*)malloc(cap), 
    .offset = 0,
    .cap = cap,
  };
};

ArenaState arenaBegin(const Arena* a) {
  return ArenaState{.restore_offset = a->offset};
}

void arenaEnd(Arena* a, ArenaState s) {
  a->offset = s.restore_offset;
}

void arenaFree(Arena *a) {
  a->cap = 0;
  a->offset = 0;
  free(a->base_ptr);
}

bool stringU8EqlZInner(const StringU8 a, const u8* bz, u64 untermd_len) {
    return (a.len == untermd_len && memcmp(a.str, bz, untermd_len) == 0); 
}

bool stringU8EqlZInner(const StringU8 a, const char* bz, u64 untermd_len) {
   return stringU8EqlZInner(a, (u8*)bz, untermd_len);
}

StringU8 stringU8FromFile(Arena* arena, FILE* in, const u64 max_read_bytes) {
  StringU8 result;
  result.cap = max_read_bytes;
  result.str = arenaAlloc(arena, u8, max_read_bytes);
  result.len = fread(result.str, 1, max_read_bytes, in);
  return result; 
}

bool isDigit(u8 ch) {
    return (ch >= '0' && ch <= '9');
}

bool isAlpha(u8 ch) {
    return ((ch >= 'a' && ch <= 'z') || 
            (ch >= 'A' && ch <= 'Z'));
}

f64 getTimeMS() {
  return ((f64)clock() / CLOCKS_PER_SEC) * 1000.0;
}
