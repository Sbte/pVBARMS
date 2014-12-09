#include "profile.h"

Profile *profiles = NULL;
int profile_amount = 0;
int profile_capacity = 0;

static Profile *find_profile(char const *name)
{
    int i;
    Profile *current_profile = NULL;
    for (i = 0; i < profile_amount; i++)
    {
        if (!strncmp(name, profiles[i].name, 31))
        {
            current_profile = &profiles[i];
            break;
        }
    }
    return current_profile;
}

void PROFILE_START_(char const *name)
{
    Profile *current_profile = find_profile(name);

    if (current_profile == NULL)
    {
        if (profile_amount == profile_capacity)
        {
            profile_capacity = (profile_capacity ? profile_capacity * 2 : 1);
            Profile *new_profiles = (Profile *)malloc(sizeof(Profile) * profile_capacity);
            memcpy(new_profiles, profiles, sizeof(Profile) * profile_amount);
            free(profiles);
            profiles = new_profiles;
        }
        current_profile = profiles + profile_amount;
        strncpy(current_profile->name, name, 31);
        parms_TimerCreate(&current_profile->timer);
        current_profile->calls = 0;
        profile_amount++;
    }
    current_profile->calls++;
    parms_TimerRestart(current_profile->timer);
}

void PROFILE_END_(char const *name)
{
    Profile *current_profile = find_profile(name);
    if (current_profile != NULL)
        parms_TimerPause(current_profile->timer);
    else
        printf("Profile %s not found\n", name);
}

void SAVE_PROFILES_(char const *filename)
{
    FILE *fp;
    fp = fopen(filename, "w");
    int i;
    Profile *profile;
    fprintf(fp, "%32s: %20s, %20s, %10s\n", "Name", "Total time", "Time per call", "Calls");
    for (i = 0; i < profile_amount; i++)
    {
        profile = &profiles[i];
        fprintf(fp, "%32s: %20f, %20f, %10d\n", profile->name, profile->timer->elapsed_time, profile->timer->elapsed_time / (double)profile->calls, profile->calls);
        parms_TimerFree(&profile->timer);
    }
    free(profiles);
    fclose(fp);
}
