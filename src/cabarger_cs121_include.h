//!
//! cabarger_cs121_include.h
//!
//! CS 121.Bolden
//! Clang 16.0.6
//! Caleb Barger
//! 09/25/23 
//! Linux x86_64
//! barg8397@vandals.uidaho.edu
//!
//! Basic types and functions for cs121 assignments this semester.
//!

#ifndef CABARGER_CS121_INCLUDE_H

#include <stdint.h>
#include <stdio.h> // FILE*

typedef uint8_t u8;
typedef int8_t i8;
typedef uint16_t u16;
typedef int16_t i16;
typedef uint32_t u32;
typedef int32_t i32;
typedef uint64_t u64;
typedef int64_t i64;

typedef float f32;
typedef double f64;

#define AssertBreak(m) (*((volatile u32*)0) = 0xCA1EB)  
#define AssertMessage(m) AssertBreak(m)
#define InvalidPath AssertMessage("invalid path")

#define KB(n) (n * 1024)
#define MB(n) (KB(n) * 1024)

#define Max(a, b) (a > b ? a : b)
#define Min(a, b) (a < b ? a : b)

struct Arena {
  u8* base_ptr;
  u64 offset; 
  u64 cap;
};

struct ArenaState {
  u64 restore_offset;
};

void arenaFree(Arena *a);
void arenaEnd(Arena* a, ArenaState s);
ArenaState arenaBegin(const Arena* a); 
Arena arenaInit(const u64 cap);
Arena subArena(Arena* arena, const u64 cap);
void* arenaAllocInner(Arena* a, const u64 requested_space); 

#define arenaAlloc(arena, T, n) \
  (T*)arenaAllocInner(arena, (u64)sizeof(T) * n) 

struct StringU8 {
    u8* str;
    u16 len; 
    u16 cap;
};

bool stringU8EqlZInner(const StringU8 a, const u8* bZ, u64 untermd_len); 
bool stringU8EqlZInner(const StringU8 a, const char* bZ, u64 untermd_len); 
#define stringU8EqlZ(a, bz) (stringU8EqlZInner(a, bz, sizeof(bz) - 1))
StringU8 stringU8FromFile(Arena* arena, FILE* in, const u64 max_read_bytes); 

bool isDigit(u8 ch); 
bool isAlpha(u8 ch); 

f64 getTimeMS();

#define CABARGER_CS121_INCLUDE_H
#endif
