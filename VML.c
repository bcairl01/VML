/// @file 	RobotLibCoreVectorMath.c
///	@brief	A set matix/vector operations geared toward robotics applications.

#include "RobotLibCore/RobotLibCoreVectorMath.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
//									  DOXYGEN GROUP DEFS										//
//////////////////////////////////////////////////////////////////////////////////////////////////
/// @defgroup rVector3
/// @defgroup rMatrix3
/// @defgroup rTransform


#if 	rDEBUG_LEVEL > 0UL
#define rCPV_CHECKSIZE(n)	if( ( (*_cpv) != n ) ) 						{ abort(); }
#define rCPM_CHECKSIZE(n,m)	if( (_cpm->nr!= n ) || (_cpm->nc!= m ) )	{ abort(); } 
#define rCPM_CHECKSQR()		if( _cpm->nr != _cpm->nc ) 					{ abort(); } 
#else
#define rCPV_CHECKSIZE(n)	
#define rCPM_CHECKSIZE(n,m)
#define rCPM_CHECKSQR()
#endif

#define VSS_1					const unsigned int*  __vss1__
#define VSS_1_PASS 				__vss1__
#define VSS_1_NROWS 			(*__vss1__)
#define VSS_1_ITR(idx)			for( idx = 0; idx < VSS_1_NROWS ; idx++ )

#define VSS_2					const unsigned int*  __vss2__
#define VSS_2_PASS 				__vss2__
#define VSS_2_NROWS 			(*__vss2__)
#define VSS_2_ITR(idx)			for( idx = 0; idx < VSS_2_NROWS ; idx++ )

#define MSS_1 					const vml_Size* __mss1__
#define MSS_1_PASS 				__mss1__
#define MSS_1_PASS_NROWS 		&(__mss1__->nr)
#define MSS_1_PASS_NCOLS 		&(__mss1__->nc)
#define MSS_1_NROWS 			__mss1__->nr
#define MSS_1_NCOLS 			__mss1__->nr
#define MSS_1_RITR(idx)			for( idx = 0; idx < MSS_1_NROWS ; idx++ )
#define MSS_1_CITR(idx)			for( idx = 0; idx < MSS_1_NCOLS ; idx++ )
#define MSS_1_SIDX(M,i,j)		MatXIDX(M,MSS_1_NROWS,i,j)
#define MSS_1_SET(var,VMLX)		var=(vml_Size*)(VMLX)
#define MSS_1_GET(var)			(const vml_Size*)(var)
#define MSS_1_GET_NROWS(var) 	((const unsigned int*)&(var->nr))
#define MSS_1_GET_NCOLS(var)	((const unsigned int*)&(var->nc))

#define MSS_2 					const vml_Size* __mss2__
#define MSS_2_PASS 				__mss2__
#define MSS_2_PASS_NROWS 		&(__mss2__->nr)
#define MSS_2_PASS_NCOLS 		&(__mss2__->nc)
#define MSS_2_NROWS 			__mss2__->nr
#define MSS_2_NCOLS 			__mss2__->nr
#define MSS_2_RITR(idx)			for( idx = 0; idx < MSS_2_NROWS ; idx++ )
#define MSS_2_CITR(idx)			for( idx = 0; idx < MSS_2_NCOLS ; idx++ )
#define MSS_2_SIDX(M,i,j)		MatXIDX(M,MSS_1_NROWS,i,j)
#define MSS_2_SET(var,VMLX)		var=(vml_Size*)(VMLX)
#define MSS_2_GET(var)			(const vml_Size*)(var)
#define MSS_2_GET_NROWS(var) 	((const unsigned int*)&(var->nr))
#define MSS_2_GET_NCOLS(var)	((const unsigned int*)&(var->nc))

//////////////////////////////////////////////////////////////////////////////////////////////////
//									      CPP-UTILS		    									//
//////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef __cplusplus

	/// Format-len
	static size_t rMTy_print_rows = 3UL;
	static size_t rMTy_print_cols = 3UL;

	/// Formatting
	std::ostream& operator| ( std::ostream& os, const rSize* fmt )
	{
		rMTy_print_rows = fmt->nr;
		rMTy_print_cols = fmt->nc;
		return os;
	}

	std::ostream& operator| ( std::ostream& os, const rU16* fmt )
	{
		rMTy_print_rows = 1.0f;
		rMTy_print_cols = *fmt;
		return os;
	}

	/// rVectorX Stream Operators
	std::ostream& operator<<( std::ostream& os, rMTy_p userp ) 
	{ 
		for( size_t idx = 0UL; idx < rMTy_print_rows; idx++ )
		{
			for( size_t jdx = 0UL; jdx < rMTy_print_cols; jdx++ )
				os << rMatXIDX(userp,rMTy_print_rows,idx,jdx) << '\t';

			if( idx!= (rMTy_print_rows-1) )
				os << '\n';
		}
		return os; 
	}

#endif



//////////////////////////////////////////////////////////////////////////////////////////////////
//									      METHODS		    									//
//////////////////////////////////////////////////////////////////////////////////////////////////


///	@ingroup Defined_Constants
/// @{

/// SIZES
const rSize _2x2 		= 	{ 2UL, 2UL };
const rSize _3x3 		= 	{ 3UL, 3UL };
const rSize _4x4 		= 	{ 4UL, 4UL };
const rSize _5x5 		= 	{ 5UL, 5UL };
const rSize _6x6 		= 	{ 6UL, 6UL };
const rSize _7x7 		= 	{ 7UL, 7UL };
const rSize _8x8 		= 	{ 8UL, 8UL };
const rSize _9x9 		= 	{ 9UL, 9UL };

const rU16 _2x1 		= 	2UL;
const rU16 _3x1 		= 	3UL;
const rU16 _4x1 		= 	4UL;
const rU16 _5x1 		= 	5UL;
const rU16 _6x1 		= 	6UL;
const rU16 _7x1 		= 	7UL;
const rU16 _8x1 		= 	8UL;
const rU16 _9x1 		= 	9UL;

const rU16 _1x2 		= 	2UL;
const rU16 _1x3 		= 	3UL;
const rU16 _1x4 		= 	4UL;
const rU16 _1x5 		= 	5UL;
const rU16 _1x6 		= 	6UL;
const rU16 _1x7 		= 	7UL;
const rU16 _1x8 		= 	8UL;
const rU16 _1x9 		= 	9UL;


/// VECTORS
rVector3 rZeroVector	= 	{ 0, 0, 0 };
rVector3 rOnesVector 	= 	{ 1, 1, 1 };
rVector3 rXhatVector	= 	{ 1, 0, 0 };
rVector3 rYhatVector	= 	{ 0, 1, 0 };
rVector3 rZhatVector	= 	{ 0, 0, 1 };
rVector3 rNOnesVector 	= 	{-1,-1,-1 };
rVector3 rNXhatVector	= 	{-1, 0, 0 };
rVector3 rNYhatVector	= 	{ 0,-1, 0 };
rVector3 rNZhatVector	= 	{ 0, 0,-1 };


/// MATRICES
rMatrix3 rIdentity3 	= 	{
								1, 0, 0,
								0, 1, 0,
								0, 0, 1,
							};
rMatrix3 rZeroMatrix3 	= 	{
								0, 0, 0,
								0, 0, 0,
								0, 0, 0,
							};

/// PLANES
rPlane	rXY_Plane		=	{ FILL3(0.0f), {0.0f, 0.0f, 1.0f} };
rPlane	rXZ_Plane		=	{ FILL3(0.0f), {0.0f, 1.0f, 0.0f} };
rPlane	rYZ_Plane		=	{ FILL3(0.0f), {1.0f, 0.0f, 0.0f} };

/// @}





///	@ingroup rVector
/// @{

	void _rVectorSetEq( rMTy_p vecOUT, rMTy_p vecIN, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = vecIN[idx];
		}
	}


	void _rVectorSetEqD2S( rMTy_p vecOUT, rMTy_dp vecIN, CPVS )
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = (rMTy)vecIN[idx];
		}
	}


	void _rVectorNegate( rMTy_p vecOUT, rMTy_p vecIN, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = -(vecIN[idx]);
		}
	
	}


	void _rVectorAdd( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, CPVS  ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = vecINA[idx] + vecINB[idx];
		}
	}


	void _rVectorAddC( rMTy_p vecOUT, rMTy_p vecINA, CPVS  ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] += vecINA[idx];
		}
	}


	void _rVectorScaleAdd( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, rF32* scale, CPVS )
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = vecINA[idx] + vecINB[idx]*(*scale);
		}		
	}


	void _rVectorScaleAddC( rMTy_p vecOUT, rMTy_p vecINA, rF32* scale, CPVS )
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] += vecINA[idx]*(*scale);
		}		
	}


	void _rVectorSub( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = vecINA[idx] - vecINB[idx];		
		}

	}


	void _rVectorSubC( rMTy_p vecOUT, rMTy_p vecINA, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] -= vecINA[idx];		
		}

	}


	void _rVectorScaleSub( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, rF32* scale, CPVS )
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = vecINA[idx] - vecINB[idx]*(*scale);
		}		
	}


	void _rVectorScaleSubC( rMTy_p vecOUT, rMTy_p vecINA, rF32* scale, CPVS )
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] -= vecINA[idx]*(*scale);
		}		
	}


	void _rVectorEMult( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = (vecINA[idx])*(vecINB[idx]);
		}

	}



	void _rVectorEDiv( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = (vecINA[idx])/(vecINB[idx]);
		}
	}



	void _rVectorESig( rMTy_p vecOUT, rMTy_p vecINA, CPVS )
	{
		unsigned int idx;
		CPV_ITR(idx)
		{
			vecOUT[idx] = 1.0f/(1.0f + rExp(-vecINA[idx]));
		}
	}



	void _rVectorScale( rMTy_p vecOUT, rMTy_p vecIN, rF32 *scale, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = (vecIN[idx])*(*scale);
		}
	}


	void _rVectorScaleC( rMTy_p vecOUT, rF32 *scale, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] *= *scale;
		}
	}


	void _rVectorScaleNR( rMTy_p vecOUT, rMTy_p vecIN, rF32 scale, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = vecIN[idx]*scale;
		}
	}


	void _rVectorScaleCNR( rMTy_p vecOUT, rF32 scale, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] *= scale;
		}
	}


	void _rVectorScaleDiv( rMTy_p vecOUT, rMTy_p vecIN, rF32 *scale, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = (vecIN[idx])/(*scale);
		}

	}



	void _rVectorShift( rMTy_p vecOUT, rMTy_p vecIN, rF32 shift, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = vecIN[idx] + shift;
		}
	}



	void _rVectorShiftC( rMTy_p vecOUT, rF32 shift, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] += shift;
		}
	}



	rF32 _rVectorDot( rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		
		rF32 dotSUM = 0;
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			dotSUM += (vecINA[idx])*(vecINB[idx]);
		}
	
		return dotSUM;
	
	}


	rF32 	_rVectorAllignment( rMTy_p vecINA, rMTy_p vecINB, CPVS )
	{
		rF32 magn = _rVectorNorm(vecINA,2.0f,CPV_PASS)*_rVectorNorm(vecINB,2.0f,CPV_PASS);

		if( magn > 0.0f )
			return (_rVectorDot(vecINA,vecINB,CPV_PASS)/magn + 1.0f)/2.0f;
		
		return -1.0f;
	}



	rF32 _rVectorNorm( rMTy_p vecIN, rF32 order, CPVS ) 
	{
		rF32 normSUM = 0;
		if(order==1) 
		{
			normSUM = rVectorDot( vecIN, vecIN );
		} 
		else 
		{
			unsigned int idx = 0;
			CPV_ITR(idx)
			{
				normSUM += rPow(vecIN[idx],order);
			}
		}
		
		return rPow(normSUM,(1/((rF32)order)));
	}


	void _rVectorProjection( rMTy_p vecOUT, rMTy_p vecIN, rPlane* plane, CPVS )
	{
		rCPV_CHECKSIZE(3UL)
		{	
			rVector3 pt;
			rVector3 dir;
			rF32 dot_out;

			rVectorSetEq(dir,plane->direction);
			rVectorSub(pt,vecIN,plane->point);
			
			dot_out = 
			rVectorDot(dir,pt);
			rVectorScaleC(dir,&dot_out);
			rVectorSub(vecOUT,vecIN,dir);
		}
	}


	rF32 _rVectorL1Diff( rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		rF32 sum = 0;
		rU08 idx = 0;

		CPV_ITR(idx)
		{
			sum += rAbs(vecINA[idx]-vecINB[idx]);
		}
		return sum;
	}



	rF32 _rVectorL2Diff( rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		rF32 sum = 0;
		rU08 idx = 0;
		CPV_ITR(idx)
			sum += rPow((vecINA[idx]-vecINB[idx]),2);

		return rSqrt(sum);
	}



	rF32 _rVectorNormalize( rMTy_p vecOUT, rMTy_p vecIN, CPVS ) 
	{
		rF32 l2inv, l2norm = rVectorNorm(vecIN,2.0f,CPV_PASS);

		if( l2norm > 1e-32 ) 
		{
			l2inv  = 1/l2norm;
		} else {
			// Epsilon Error
			return 0;
		}

		rVectorScale(vecOUT,vecIN,&l2inv,CPV_PASS);

		return(l2norm);
	}



	rF32 _rVectorNormalizeC( rMTy_p vecOUT, CPVS )
	{
		rF32 l2inv, l2norm = rVectorNorm(vecOUT,2.0f,CPV_PASS);

		if( l2norm > M_WORKING_PRECISION ) 
		{
			l2inv  = 1/l2norm;
		} else {
			// Epsilon Error
			return 0;
		}

		rVectorScale(vecOUT,vecOUT,&l2inv,CPV_PASS);

		return(l2norm);
	}



	rF32 _rVectorNormalizeDiff( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		
		rF32 l2inv, l2norm = rVectorL2Diff(vecINA,vecINB,CPV_PASS);

		if( l2norm > 1e-32 ) 
		{
			l2inv = 1/l2norm;
		} else {
			// Epsilon Error
			return 0;
		}

		rVectorSub(vecOUT,vecINA,vecINB,CPV_PASS);
		
		rVectorScale(vecOUT,vecOUT,&l2inv,CPV_PASS);

		return l2norm;
	}



	void _rVectorMaskMax( rMTy_p vecOUT, rMTy_p vecIN, CPVS )
	{
		unsigned int idx;
		unsigned int mdx;
		for( idx = 1UL ; idx < *CPV_PASS ; idx++ )
		{
			if( vecIN[idx] > vecIN[idx-1UL])
				mdx = idx;
			else
				mdx = idx-1;
		}
		CPV_ITR(idx)
		{
			if(idx!=mdx)
				vecOUT[idx] = 0.0f;
			else
				vecOUT[idx] = vecIN[idx];
		}
	}


	void _rVectorMaskMaxC( rMTy_p vecOUT, CPVS )
	{
		unsigned int idx;
		unsigned int mdx;
		for( idx = 1UL ; idx < *CPV_PASS ; idx++ )
		{
			if( vecOUT[idx] > vecOUT[idx-1UL])
				mdx = idx;
			else
				mdx = idx-1;
		}
		CPV_ITR(idx)
		{
			if(idx!=mdx)
				vecOUT[idx] = 0.0f;
		}
	}



	void _rVectorMaskMin( rMTy_p vecOUT, rMTy_p vecIN, CPVS )
	{
		unsigned int idx;
		unsigned int mdx;
		for( idx = 1UL ; idx < *CPV_PASS ; idx++ )
		{
			if( vecIN[idx] < vecIN[idx-1UL])
				mdx = idx;
			else
				mdx = idx-1;
		}
		CPV_ITR(idx)
		{
			if(idx!=mdx)
				vecOUT[idx] = 0.0f;
			else
				vecOUT[idx] = vecIN[idx];
		}
	}


	void _rVectorMaskMinC( rMTy_p vecOUT, CPVS )
	{
		unsigned int idx;
		unsigned int mdx;
		for( idx = 1UL ; idx < *CPV_PASS ; idx++ )
		{
			if( vecOUT[idx] < vecOUT[idx-1UL])
				mdx = idx;
			else
				mdx = idx-1;
		}
		CPV_ITR(idx)
		{
			if(idx!=mdx)
				vecOUT[idx] = 0.0f;
		}
	}


	void _rVectorCross( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{		
			rVector3 tempA;
			rVector3 tempB;

			rVectorSetEq(tempA,vecINA);
			rVectorSetEq(tempB,vecINB);

			vecOUT[0] =  (tempA[1])*(tempB[2]) - (tempA[2])*(tempB[1]);
			vecOUT[1] =  (tempA[2])*(tempB[0]) - (tempA[0])*(tempB[2]);
			vecOUT[2] =  (tempA[0])*(tempB[1]) - (tempA[1])*(tempB[0]);
		}
	}


	rF32 _rVectorTriple( rMTy_p vecINA, rMTy_p vecINB, rMTy_p vecINC, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{	
			return
				((vecINA[1])*(vecINB[2]) - (vecINA[2])*(vecINB[1]))*vecINC[0] + 
				((vecINA[2])*(vecINB[0]) - (vecINA[0])*(vecINB[2]))*vecINC[1] + 
				((vecINA[0])*(vecINB[1]) - (vecINA[1])*(vecINB[0]))*vecINC[2];
		}
		return 0.0f;
	}



	void _rVectorAverage( rMTy_p vecOUT, rMTy_p vecINA, rMTy_p vecINB, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx) 
		{
			vecOUT[idx] = (vecINA[idx]+vecINB[idx])/2;	
		}
	}


	void _rVectorAbs( rMTy_p vecOUT, rMTy_p vecINA, CPVS) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx) 
		{
			vecOUT[idx] = rAbs(vecINA[idx]);	
		}
	}


	void _rVectorAbsC( rMTy_p vecOUT, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx) 
		{
			vecOUT[idx] = rAbs(vecOUT[idx]);	
		}
	}


	void _rVectorXRotation( rMTy_p vecOUT, rMTy_p vecIN, rF32 *angle, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{		
			rVector3 tempA;
			rVectorSetEq(tempA,vecIN);
			vecOUT[0] = tempA[0];
			vecOUT[1] = rCos(*angle)*tempA[1] - rSin(*angle)*tempA[2];
			vecOUT[2] = rSin(*angle)*tempA[1] + rCos(*angle)*tempA[2];
		}
	}




	void _rVectorYRotation( rMTy_p vecOUT, rMTy_p vecIN, rF32 *angle, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{				
			rVector3 tempA;
			rVectorSetEq(tempA, vecIN);
			
			vecOUT[0] = rSin(*angle)*tempA[2] + rCos(*angle)*tempA[0];
			vecOUT[1] = tempA[1];
			vecOUT[2] = rCos(*angle)*tempA[2] - rSin(*angle)*tempA[0];
		}
	}




	void _rVectorZRotation( rMTy_p vecOUT, rMTy_p vecIN, rF32 *angle, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{				
			rVector3 tempA;
			rVectorSetEq(tempA, vecIN);
			
			vecOUT[0] = rCos(*angle)*tempA[0] - rSin(*angle)*tempA[1];
			vecOUT[1] = rSin(*angle)*tempA[0] + rCos(*angle)*tempA[1];
			vecOUT[2] = tempA[2]; 
		}
	}




	void _rVectorXRotationC( rMTy_p vecOUT, rF32 *angle, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{			
			rVector3 tempA;
			rVectorSetEq(tempA,vecOUT);
			
			vecOUT[1] = rCos(*angle)*tempA[1] - rSin(*angle)*tempA[2];
			vecOUT[2] = rSin(*angle)*tempA[1] + rCos(*angle)*tempA[2];
		}
	}




	void _rVectorYRotationC( rMTy_p vecOUT, rF32 *angle, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{	
			rVector3 tempA;
			rVectorSetEq(tempA, vecOUT);
			
			vecOUT[0] = rSin(*angle)*tempA[2] + rCos(*angle)*tempA[0];
			vecOUT[2] = rCos(*angle)*tempA[2] - rSin(*angle)*tempA[0];
		}
	}




	void _rVectorZRotationC( rMTy_p vecOUT, rF32 *angle, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{	
			rVector3 tempA;
			rVectorSetEq(tempA, vecOUT, CPV_PASS);
			
			vecOUT[0] = rCos(*angle)*tempA[0] - rSin(*angle)*tempA[1];
			vecOUT[1] = rSin(*angle)*tempA[0] + rCos(*angle)*tempA[1];
		}
	}


	void _rVectorLoadPolar( rMTy_p vecOUT, rF32* theta, rF32* alpha, CPVS )
	{
		rCPV_CHECKSIZE(3UL)
		{	
			rVec_X(vecOUT) =-rSin(*theta)*rCos(*alpha);
			rVec_Y(vecOUT) = rCos(*theta)*rCos(*alpha);
			rVec_Z(vecOUT) = rSin(*alpha);
		}
	}


	void 	_rVectorLoadPolar90a( rMTy_p vecOUT, rF32* theta, rF32* alpha, CPVS )
	{
		rCPV_CHECKSIZE(3UL)
		{	
			rVec_X(vecOUT) =-rSin(*theta)*rCos(*alpha + M_PI1_2);
			rVec_Y(vecOUT) = rCos(*theta)*rCos(*alpha + M_PI1_2);
			rVec_Z(vecOUT) = rSin(*alpha + M_PI1_2);
		}
	}


	void 	_rVectorLoadPolar90t( rMTy_p vecOUT, rF32* theta, rF32* alpha, CPVS )
	{
		rCPV_CHECKSIZE(3UL)
		{	
			rVec_X(vecOUT) =-rSin(*theta + M_PI1_2)*rCos(*alpha );
			rVec_Y(vecOUT) = rCos(*theta + M_PI1_2)*rCos(*alpha);
			rVec_Z(vecOUT) = rSin(*alpha);
		}
	}


	void _rVectorLoadZeros( rMTy_p vecOUT, CPVS ) 
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = 0;
		}
	}


	void 	_rVectorLoadRand( rMTy_p vecOUT, CPVS )
	{
		unsigned int idx = 0;
		CPV_ITR(idx)
		{
			vecOUT[idx] = ((rMTy)(rand()%10000UL)/1e5f);
		}
	}


	void _rVectorLoadXAxis( rMTy_p vecOUT, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{
			vecOUT[0] = 1.0f;
			vecOUT[1] = 0.0f;
			vecOUT[2] = 0.0f;
		}
	}


	void _rVectorLoadYAxis( rMTy_p vecOUT, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{
			vecOUT[0] = 0.0f;
			vecOUT[1] = 1.0f;
			vecOUT[2] = 0.0f;
		}
	}


	void _rVectorLoadZAxis( rVector3 vecOUT, CPVS ) 
	{
		rCPV_CHECKSIZE(3UL)
		{
			vecOUT[0] = 0.0f;
			vecOUT[1] = 0.0f;
			vecOUT[2] = 1.0f;
		}
	}


	void 	_rVectorLoadAll( rMTy_p vecOUT, rF32* num, CPVS )
	{
		unsigned int idx;
		CPV_ITR(idx)
		{
			vecOUT[idx]= *num;
		}
	}


	void rVector3RPY( rVector3 vecOUT, rMatrix3 rotation )
	{
		/// @note Refer to http://planning.cs.uiuc.edu/node103.html
		rVec_Z(vecOUT) = rArcTan2( rMat3_21(rotation), rMat3_11(rotation));
		rVec_Y(vecOUT) = rArcTan2(-rMat3_31(rotation), rSqrt(rPow2(rMat3_32(rotation) + rMat3_33(rotation))) );
		rVec_X(vecOUT) = rArcTan2( rMat3_32(rotation), rMat3_33(rotation));
	}

/// @}










/// @ingroup rMatrix3
/// @{



	void _rMatrixSetEq( rMatrix3 matOUT, rMatrix3 matIN, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) = CPM_SIDX(matIN,idx,jdx);
			}
		}
	}



	void _rMatrixScale( rMTy_p matOUT, rMTy_p matIN, rF32 *scale, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) = CPM_SIDX(matIN,idx,jdx)*(*scale);
			}
		}
	}



	void _rMatrixTranspose( rMTy_p matOUT, rMTy_p matIN, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,jdx,idx) = CPM_SIDX(matIN,idx,jdx);
			}
		}
	}



	void _rMatrixInverse( rMTy_p matOUT, rMTy_p matIN, CPMS ) 
	{
		rCPM_CHECKSQR()
		{

		}
		/// @todo
	}



	void _rMatrixDeterminant( rMTy_p matOUT, rMTy_p matIN, CPMS ) 
	{
		rCPM_CHECKSQR()
		{


		}
		/// @todo
	}



	void _rMVS_Mult( rMTy_p vecOUT, rMTy_p matIN, rMTy_p vecIN, CPMS ) 
	{
		rTHROW_ON_MATCH(vecOUT,vecIN);

		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			vecOUT[idx] = 0;
			CPM_CITR(jdx)
			{
				vecOUT[idx] += CPM_SIDX(matIN,idx,jdx)*(vecIN[jdx]);
			}
		}
	}


	void 	_rVMS_Mult( rMTy_p vecOUT, rMTy_p vecIN, rMTy_p matIN, CPMS )
	{
		rTHROW_ON_MATCH(vecOUT,vecIN);

		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			vecOUT[idx] = 0;
			CPM_CITR(jdx)
			{
				vecOUT[idx] += (vecIN[jdx])*CPM_SIDX(matIN,jdx,idx);
			}
		}
	}



	void _rMMS_Mult(  rMTy_p matOUT, rMTy_p matINA, rMTy_p matINB, CPMS ) 
	{
		rTHROW_ON_MATCH(matOUT,matINA);
		rTHROW_ON_MATCH(matOUT,matINB);

		unsigned int idx, jdx, kdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,jdx,idx) = 0;
				CPM_RITR(kdx)
				{
					CPM_SIDX(matOUT,jdx,idx) += CPM_SIDX(matINA,jdx,kdx)*CPM_SIDX(matINB,kdx,idx);
				}
			}
		}
	}



	void _rMatrixLoadIdentity( rMTy_p matOUT, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				if(idx==jdx){
					CPM_SIDX(matOUT,idx,jdx) = 1;
				} else {
					CPM_SIDX(matOUT,idx,jdx) = 0;
				}
			}
		}
	}




	void _rMatrixLoadZeros( rMTy_p matOUT, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) = 0;
			}
		}
	}




	void _rMatrixLoadOnes( rMTy_p matOUT, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) = 1.0f;
			}
		}
	}


	
	RLCORE_API void		_rMatrixLoadRand( rMTy_p matOUT, CPMS )
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) = ((rMTy)(rand()%10000UL)/1e5f);;
			}
		}
	}



	void 	_rMatrixAdd( rMTy_p matOUT, rMTy_p matINA, rMTy_p matINB, CPMS )
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) = CPM_SIDX(matINA,idx,jdx) + CPM_SIDX(matINB,idx,jdx);
			}
		}
	}



	void 	_rMatrixSub( rMTy_p matOUT, rMTy_p matINA, rMTy_p matINB, CPMS )
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) = CPM_SIDX(matINA,idx,jdx) - CPM_SIDX(matINB,idx,jdx);
			}
		}
	}



	void _rMatrixAddC( rMTy_p matOUT, rMTy_p matIN, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) += CPM_SIDX(matIN,idx,jdx);
			}
		}
	}




	void _rMatrixSubC( rMTy_p matOUT, rMTy_p matIN, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) -= CPM_SIDX(matIN,idx,jdx);
			}
		}
	}



	void _rMatrixScaleAddC( rMTy_p matOUT, rMTy_p matIN, rF32* scale, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) += CPM_SIDX(matIN,idx,jdx)*(*scale);
			}
		}
	}




	void _rMatrixScaleSubC( rMTy_p matOUT, rMTy_p matIN, rF32* scale, CPMS ) 
	{
		unsigned int idx, jdx;
		CPM_RITR(idx)
		{
			CPM_CITR(jdx)
			{
				CPM_SIDX(matOUT,idx,jdx) -= CPM_SIDX(matIN,idx,jdx)*(*scale);
			}
		}
	}




	void rMatrix3LoadSkew( rMatrix3 matOUT, rVector3 vecIN ) 
	{

		// Zero Diagonal
		unsigned int idx = 0;
		for( ; idx < 3; idx++ ) 
		{
			rMat3IDX(matOUT,idx,idx) = 0;
		}

		// Upper Triangle
		rMat3IDX(matOUT,0,1) = -(vecIN[0]);
		rMat3IDX(matOUT,0,2) =  (vecIN[1]);
		rMat3IDX(matOUT,1,2) = -(vecIN[2]);

		// Lower Triangle
		rMat3IDX(matOUT,1,0) = -rMat3IDX(matOUT,0,1);
		rMat3IDX(matOUT,2,0) =  rMat3IDX(matOUT,0,2);
		rMat3IDX(matOUT,2,1) = -rMat3IDX(matOUT,1,2);	
	}


	void rMatrix3LoadZRot( rMatrix3 matOUT, rF32 angle ) 
	{
		rMatrixLoadIdentity(matOUT);
		rMat3IDX(matOUT,0,0) = rCos(angle);
		rMat3IDX(matOUT,0,1) =-rSin(angle);
		rMat3IDX(matOUT,1,0) = rSin(angle);
		rMat3IDX(matOUT,1,1) = rCos(angle);
	}


	void rMatrix3LoadYRot( rMatrix3 matOUT, rF32 angle ) 
	{
		rMatrixLoadIdentity(matOUT);
		rMat3IDX(matOUT,0,0) = rCos(angle);
		rMat3IDX(matOUT,0,2) = rSin(angle);
		rMat3IDX(matOUT,2,0) =-rSin(angle);
		rMat3IDX(matOUT,2,2) = rCos(angle);
	}



	void rMatrix3LoadXRot( rMatrix3 matOUT, rF32 angle ) 
	{
		rMatrixLoadIdentity(matOUT);
		rMat3IDX(matOUT,1,1) = rCos(angle);
		rMat3IDX(matOUT,1,2) =-rSin(angle);
		rMat3IDX(matOUT,2,1) = rSin(angle);
		rMat3IDX(matOUT,2,2) = rCos(angle);
	}


	void rMatrix3LoadRPY( rMatrix3 matOUT, rVector3 vecRPY ) 
	{
		// R1
		rMat3IDX(matOUT,0,0) = rCos(rVec_Y(vecRPY))*rCos(rVec_Z(vecRPY));
		rMat3IDX(matOUT,0,1) = rCos(rVec_Z(vecRPY))*rSin(rVec_X(vecRPY)) * rSin(rVec_Y(vecRPY)) - rCos(rVec_X(vecRPY))*rSin(rVec_Z(vecRPY));
		rMat3IDX(matOUT,0,2) = rSin(rVec_X(vecRPY))*rSin(rVec_Z(vecRPY)) + rCos(rVec_X(vecRPY)) * rCos(rVec_Z(vecRPY))*rSin(rVec_Y(vecRPY));
		// R2
		rMat3IDX(matOUT,1,0) = rCos(rVec_Y(vecRPY))*rSin(rVec_Z(vecRPY));
		rMat3IDX(matOUT,1,1) = rCos(rVec_X(vecRPY))*rCos(rVec_Z(vecRPY)) + rSin(rVec_X(vecRPY)) * rSin(rVec_Y(vecRPY))*rSin(rVec_Z(vecRPY));
		rMat3IDX(matOUT,1,2) = rCos(rVec_X(vecRPY))*rSin(rVec_Y(vecRPY)) * rSin(rVec_Z(vecRPY)) - rCos(rVec_Z(vecRPY))*rSin(rVec_X(vecRPY));
		// R3
		rMat3IDX(matOUT,2,0) =-rSin(rVec_Y(vecRPY));
		rMat3IDX(matOUT,2,1) = rCos(rVec_Y(vecRPY))*rSin(rVec_X(vecRPY));
		rMat3IDX(matOUT,2,2) = rCos(rVec_X(vecRPY))*rCos(rVec_Y(vecRPY));
	}

/// @}//




///	@ingroup rIPInternals
/// @{

	static rMTy_p 	rIP_mem 	= NULL;
	static rSize 	rIP_size    = {0UL,0UL};

	void _rIPResolveMatrix( CPMS )
	{
		if(rIP_mem==NULL)
		{
			rIP_mem = NEW_MATRIX(*CPM_PASS);
			rIP_size.nr = CPM_PASS->nr;
			rIP_size.nc = CPM_PASS->nc;
		}
	}


	rBool 	_rIPInUse()
	{
		return rIP_mem!=NULL;
	}


	rMTy_p	_rIPGet()
	{
		return rIP_mem;
	}


	rMTy_p	_rIPClear()
	{
		free(rIP_mem);
		rIP_mem = NULL;
		return rIP_mem;
	}

/// @}



///	@ingroup rIPMatrix
/// @{


	rMTy_p 	_rIPMatrixTranspose( rMTy_p matIN, CPMS )
	{
		_rIPResolveMatrix(CPM_PASS);
		_rMatrixTranspose(rIP_mem,matIN,CPM_PASS);
		return rIP_mem;
	}



/// @}






///	@ingroup rTransform
/// @{


	void rTransformSetEq( rTransform *homOUT, rTransform *homIN ) 
	{
		rMatrixSetEq(homOUT->rotation,homIN->rotation);
		rVectorSetEq(homOUT->translation,homIN->translation);
	}


	void rTransformLoadIdentity( rTransform *homOUT ) 
	{
		rMatrixLoadIdentity( homOUT->rotation);
		rVectorLoadZeros( homOUT->translation);
	}
	

	void rTransformLoadDHParam( rTransform *homOUT, rF32 xoff, rF32 alpha, rF32 theta, rF32 zoff ) 
	{
		rMatrixLoadIdentity( homOUT->rotation );
		
		///ROW-1
		rMat3IDX(homOUT->rotation,0,0)	= rCos(theta);
		rMat3IDX(homOUT->rotation,0,1)	=-rSin(theta)*rCos(alpha);
		rMat3IDX(homOUT->rotation,0,2)	= rSin(theta)*rSin(alpha);
		homOUT->translation[0]			= rCos(theta)*xoff;
		
		///ROW-2
		rMat3IDX(homOUT->rotation,1,0)	= rSin(theta);
		rMat3IDX(homOUT->rotation,1,1)	= rCos(theta)*rCos(alpha);
		rMat3IDX(homOUT->rotation,1,2)	=-rCos(theta)*rSin(alpha);
		homOUT->translation[1]			= rSin(theta)*xoff;

		///ROW-3
		rMat3IDX(homOUT->rotation,2,1)	= rSin(alpha);
		rMat3IDX(homOUT->rotation,2,2)	= rCos(alpha);
		homOUT->translation[2]			= zoff;
	}


	void rHV_Mult( rVector3 vecOUT, rTransform *homIN, rVector3 vecIN ) 
	{
		
		// Rotation
		_rMVS_Mult( vecOUT, homIN->rotation, vecIN, CP3x3);

		// translation
		rVectorAdd( vecOUT, homIN->translation, vecIN, CP3x1);

	}


	void rIHV_Mult( rVector3 vecOUT, rTransform *homIN, rVector3 vecIN  ) 
	{
		rMatrix3 R_inverse, safe_copy;

		// Invert Rotation rMatrix
		rMatrixTranspose( R_inverse, homIN->rotation, CP3x3); 

		// translation (negative)
		rVectorSub( vecOUT, vecIN, homIN->translation, CP3x1);

		// Rotation By Inverse
		rMatrixSetEq(safe_copy,vecOUT);
		_rMVS_Mult( vecOUT, R_inverse, safe_copy, CP3x3);
	}


	void rHH_Mult( rTransform *homOUT, rTransform *homINA, rTransform *homINB ) 
	{

		// Rotation
		_rMMS_Mult( homOUT->rotation,    homINA->rotation, homINB->rotation, CP3x3 );
		_rMVS_Mult( homOUT->translation, homINA->rotation, homINB->translation, CP3x3);

		// translation
		rVectorAdd( homOUT->translation, homOUT->translation, homINA->translation, CP3x1 );

	}


	void rIHH_Mult( rTransform *homOUT, rTransform *homINA, rTransform *homINB  ) 
	{

		rMatrix3 RA_inverse;

		// Invert Rotation rMatrix
		rMatrixTranspose( RA_inverse, homINA->rotation ); 

		// translation (negative)
		rVectorSub( homOUT->translation, homINB->translation, homINA->translation );

		// Rotation By Inverse Rotation-A
		_rMVS_Mult( homOUT->translation, RA_inverse, homOUT->translation, CP3x3);
		_rMMS_Mult( homOUT->rotation,    RA_inverse, homINB->rotation,    CP3x3);
	}


/// @}



#ifdef rUSING_OPENGL
void rVector3GLRenderAsVertex( rVector3 vec )
{
	glVertex3f(vec[0UL],vec[1UL],vec[2UL]);
}


void rVector3GLRenderAsNormal( rVector3 vec )
{
	glNormal3f(vec[0UL],vec[1UL],vec[2UL]);
}
#endif




#ifdef rUSING_SIMD
void rMatrix3_to_SIMD( rMatrix3 Mat, double simd[12] ) 
{
	for( size_t idx = 0; idx < 12; idx++ )
		simd[idx] = 0;
	for( size_t idx = 0; idx<3; idx++ ) 
	{
		for( size_t jdx = 0; jdx<3; jdx++ ) 
		{
			simd[idx*4+jdx] = rMat3IDX(Mat,idx,jdx);
		}
	}
}


void SIMD_to_rMatrix3( double simd[12], rMatrix3 Mat ) 
{
	for( size_t idx = 0; idx<3; idx++ ) 
	{
		for( size_t jdx = 0; jdx<3; jdx++ ) 
		{
			rMat3IDX(Mat,idx,jdx) = simd[idx*4+jdx];
		}
	}
}


void SIMDLoadZRot( double simd[12], double a   ) 
{
	simd[0] = +rCos(a);	simd[1] = -rSin(a);		simd[2] = 0.0;		simd[3] = 0;
	simd[4] = +rSin(a);	simd[5] = +rCos(a);		simd[6] = 0.0;		simd[7] = 0;
	simd[8] = 0.0;		simd[9] = 0.0;;			simd[10]= 1.0;		simd[11]= 0;
}


void SIMDLoadYRot( double simd[12], double a   ) 
{
	simd[0] = +rCos(a);	simd[1] = 0.0;			simd[2] = +rSin(a);	simd[3] = 0;
	simd[4] = 0.0;		simd[5] = 1.0;			simd[6] = 0.0;		simd[7] = 0;
	simd[8] = -rSin(a);	simd[9] = 0.0;;			simd[10]= +rCos(a);	simd[11]= 0;
}


void SIMDLoadXRot( double simd[12], double a   ) 
{
	simd[0] = +1.0;		simd[1] = 0.0;			simd[2] = 0.0;		simd[3] = 0;
	simd[4] = 0.0;		simd[5] = +rCos(a);		simd[6] = -rSin(a);	simd[7] = 0;
	simd[8] = 0.0;		simd[9] = +rSin(a);		simd[10]= +rCos(a);	simd[11]= 0;
}


void SIMDLoadIdentity( double simd[12] ) 
{
	simd[0] = 1.0;		simd[1] = 0.0;			simd[2] = 0.0;		simd[3] = 0;
	simd[4] = 0.0;		simd[5] = 1.0;			simd[6] = 0.0;		simd[7] = 0;
	simd[8] = 0.0;		simd[9] = 0.0;			simd[10]= 1.0;		simd[11]= 0;
}


const double SIMD_90_X[12] = 
{	+1.0,		+0.0,		+0.0,		+0.0,
	+0.0,		+0.0,		-1.0,		+0.0,
	+0.0,		+1.0,		+0.0,		+0.0,
};

const double SIMD_90_Y[12] = 
{	+0.0,		+0.0,		+1.0,		+0.0,
	+0.0,		+1.0,		-0.0,		+0.0,
	-1.0,		+0.0,		+0.0,		+0.0,
};

const double SIMD_90_Z[12] = 
{	+0.0,		-1.0,		+0.0,		+0.0,
	+1.0,		+0.0,		-0.0,		+0.0,
	+0.0,		+0.0,		+1.0,		+0.0,
};


const double SIMD_N90_X[12] = 
{	+1.0,		+0.0,		+0.0,		+0.0,
	+0.0,		+0.0,		+1.0,		+0.0,
	+0.0,		-1.0,		+0.0,		+0.0,
};

const double SIMD_N90_Y[12] = 
{	+0.0,		+0.0,		-1.0,		+0.0,
	+0.0,		+1.0,		-0.0,		+0.0,
	+1.0,		+0.0,		+0.0,		+0.0,
};

const double SIMD_N90_Z[12] = 
{	+0.0,		+1.0,		+0.0,		+0.0,
	-1.0,		+0.0,		-0.0,		+0.0,
	+0.0,		+0.0,		+1.0,		+0.0,
};
#endif

