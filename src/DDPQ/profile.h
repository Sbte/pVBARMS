#ifndef __ITSOL_INCLUDED_PROFILE_H__
#define __ITSOL_INCLUDED_PROFILE_H__

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../include/parms_timer_impl.h"

typedef struct
{
    char name[32];
    int calls;
    parms_Timer timer;
} Profile;

extern Profile *profiles;
extern int profile_amount;
extern int profile_capacity;

void PROFILE_START_(char const *name);
void PROFILE_END_(char const *name);
void SAVE_PROFILES_(char const *filename);

#ifndef PROFILE_ON

#define PROFILE_START(a)
#define PROFILE_END(a)
#define SAVE_PROFILES(a)

#else

#define PROFILE_START(a) PROFILE_START_(a)
#define PROFILE_END(a) PROFILE_END_(a)
#define SAVE_PROFILES(a) SAVE_PROFILES_(a)

#endif

#endif
