#ifndef _included_timer_h
#define _included_timer_h

#include "Statistics.h"

#define STATS_FIRSTPASS	(0)
#define STATS_REDUCTION (1)
#define STATS_COARSEGRID (2)
#define STATS_SETBOUNDARY	(3)
#define STATS_FINALPASS (4)
#define STATS_TOTAL (5)

#define NUM_STATS	(6)

const char *StatsNames[NUM_STATS] = 
{
  "First Pass",
  "Reduction",
  "Coarse Grid",
  "Set Boundary",
  "Final Pass",
  "Total"
};


#define STATS_START(name)	StatMaster.Start(name)
#define STATS_STOP(name)	StatMaster.Stop(name)
#define STATS_RESET()		StatMaster.Reset()
#define STATS_REPORT(name)	StatMaster.Report(name)

#endif
