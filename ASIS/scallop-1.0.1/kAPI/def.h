#ifndef _def_h
#define _def_h

/*

aliases for various dimension-dependent types

*/


#define POINT 		PointNDIM
#define REGION 		RegionNDIM
#define FLOORPLAN 	FloorPlanNDIM
#define GHOSTPLAN 	GhostPlanNDIM
#define GHOSTITERATOR 	GhostIteratorNDIM
#define S_GRID		s_gridNDIM
#define P_GRID		p_gridNDIM
#define GRID 		GridNDIM<double>
#define XARRAY 		XArrayNDIM<s_gridNDIM>
#define MOTIONPLAN	MotionPlanNDIM
#define MOVER		VectorMoverNDIM<s_gridNDIM,double>
#define ADDER		adderNDIM
#define SUBTRACTOR	subtractorNDIM
#define P_DOMAIN	p_domainNDIM

#endif
