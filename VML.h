/// @file 	VML.h
///	@brief	Vector-Math Library. Meant to be like eigen for C.
///
///	@todo
///	
/// @defgroup vml_Vector
/// @defgroup vml_Matrix
/// @defgroup vml_Transform


#ifndef VECTORMATHLIB_H
#define VECTORMATHLIB_H
#define VML_API


#ifndef VML_DEBUG_LEVEL
#define VML_DEBUG_LEVEL 0UL
#endif 


#ifndef VML_WITH_CONSTANTS
#define VML_WITH_CONSTANTS  1
#endif


#ifndef	VML_WITH_SIZEDEFS
#define VML_WITH_SIZEDEFS	1
#endif


#ifndef VML_WITH_OPENGL
#define VML_WITH_OPENGL  	0
#endif


#if 	VML_DEBUG_LEVEL > 0UL
#define VML_THROW_ON_MATCH(var1,var2) 	if(var1==var2) { abort(); }
#else
#define VML_THROW_ON_MATCH(var1,var2)
#endif


typedef float	vmlTy;
typedef double	vmlTy_d;
typedef float*	vmlTy_p;
typedef double*	vmlTy_dp;


#if __cplusplus
	#include 	<cmath>
	#include 	<iostream>
	#include 	<iomanip>
	using namespace std;
#else
	#include 	<math.h>
#endif

#if	VML_WITH_OPENGL
	#if _WIN32 || _WIND64
	#include	"Windows.h"
	#endif
	#include	<GL/GL.h>
	#include	<GL/freeglut.h>
#endif

#include		<stdlib.h>




/// @ingroup vml_Size
/// @{

	typedef struct xvml_Size
	{
		unsigned int nr;
		unsigned int nc;
	} vml_Size;


#if	VML_WITH_SIZEDEFS
	extern const vml_Size 			_2x2;
	extern const vml_Size 			_3x3;
	extern const vml_Size 			_4x4;
	extern const vml_Size 			_5x5;
	extern const vml_Size 			_6x6;
	extern const vml_Size 			_7x7;
	extern const vml_Size 			_8x8;
	extern const vml_Size 			_9x9;

	extern const unsigned int 		_2x1;
	extern const unsigned int 		_3x1;
	extern const unsigned int 		_4x1;
	extern const unsigned int 		_5x1;
	extern const unsigned int 		_6x1;
	extern const unsigned int 		_7x1;
	extern const unsigned int 		_8x1;
	extern const unsigned int 		_9x1;

	extern const unsigned int 		_1x2;
	extern const unsigned int 		_1x3;
	extern const unsigned int 		_1x4;
	extern const unsigned int 		_1x5;
	extern const unsigned int 		_1x6;
	extern const unsigned int 		_1x7;
	extern const unsigned int 		_1x8;
	extern const unsigned int 		_1x9;


	#define VML2x2					(const vml_Size*)(&_2x2)
	#define VML3x3					(const vml_Size*)(&_3x3)
	#define VML4x4					(const vml_Size*)(&_4x4)
	#define VML5x5					(const vml_Size*)(&_5x5)
	#define VML6x6					(const vml_Size*)(&_6x6)
	#define VML7x7					(const vml_Size*)(&_7x7)
	#define VML8x8					(const vml_Size*)(&_8x8)
	#define VML9x9					(const vml_Size*)(&_9x9)

	#define VML2x1					(const unsigned int*)(&_2x1)
	#define VML3x1					(const unsigned int*)(&_3x1)
	#define VML4x1					(const unsigned int*)(&_4x1)
	#define VML5x1					(const unsigned int*)(&_5x1)
	#define VML6x1					(const unsigned int*)(&_6x1)
	#define VML7x1					(const unsigned int*)(&_7x1)
	#define VML8x1					(const unsigned int*)(&_8x1)
	#define VML9x1					(const unsigned int*)(&_9x1)

	#define VML1x2					(const unsigned int*)(&_1x2)
	#define VML1x3					(const unsigned int*)(&_1x3)
	#define VML1x4					(const unsigned int*)(&_1x4)
	#define VML1x5					(const unsigned int*)(&_1x5)
	#define VML1x6					(const unsigned int*)(&_1x6)
	#define VML1x7					(const unsigned int*)(&_1x7)
	#define VML1x8					(const unsigned int*)(&_1x8)
	#define VML1x9					(const unsigned int*)(&_1x9)
#endif


#define VSS_1						const unsigned int*  __vss1__
#define VSS_2						const unsigned int*  __vss2__
#define MSS_1 						const vml_Size* __mss1__
#define MSS_2 						const vml_Size* __mss2__


/// @}




///	@ingroup vml_VectorN
/// @{
	typedef vmlTy_p vml_Vector;
 	typedef vmlTy 	vml_Vector2[2UL];  		///< A 3x1 vml_Vector type
 	typedef vmlTy 	vml_Vector3[3UL];  		///< A 3x1 vml_Vector type
  	typedef vmlTy 	vml_Vector4[4UL];  		///< A 4x1 vml_Vector type
   	typedef vmlTy 	vml_Vector5[5UL];  		///< A 5x1 vml_Vector type
   	typedef vmlTy 	vml_Vector6[6UL];  		///< A 6x1 vml_Vector type
   	typedef vmlTy 	vml_Vector7[7UL];  		///< A 7x1 vml_Vector type
   	typedef vmlTy 	vml_Vector8[8UL];  		///< A 8x1 vml_Vector type
   	typedef vmlTy 	vml_Vector9[9UL];  		///< A 9x1 vml_Vector type

   	#define NEW_COLVECTOR(size)				(vmlTy_p)(malloc(sizeof(vmlTy)*(size).nr))
   	#define NEW_ROWVECTOR(size)				(vmlTy_p)(malloc(sizeof(vmlTy)*(size).nc))
   	#define FREE_VECTOR(ptr)				free((vmlTy_p)ptr)
/// @}



///	@ingroup vml_MatrixN
/// @{
	typedef vmlTy_p vml_Matrix;
 	typedef vmlTy 	vml_Matrix2[4UL];		///< A 2x2 vml_Matrix type
 	typedef vmlTy 	vml_Matrix3[9UL];		///< A 3x3 vml_Matrix type
 	typedef vmlTy 	vml_Matrix4[16UL];		///< A 4x4 vml_Matrix type
 	typedef vmlTy 	vml_Matrix5[25UL];		///< A 5x5 vml_Matrix type
 	typedef vmlTy 	vml_Matrix6[36UL];		///< A 6x6 vml_Matrix type
 	typedef vmlTy 	vml_Matrix7[49UL];		///< A 7x7 vml_Matrix type
 	typedef vmlTy 	vml_Matrix8[64UL];		///< A 8x8 vml_Matrix type
 	typedef vmlTy 	vml_Matrix9[81UL];		///< A 9x9 vml_Matrix type

   	#define NEW_MATRIX(size)				(vmlTy_p)(malloc(sizeof(vmlTy)*(size).nr*(size).nc))
   	#define FREE_MATRIX(ptr)				free((vmlTy_p)ptr)
/// @}



///	@ingroup vml_Transform
/// @{
 	/// @brief	4x4 homogeneous vml_Transform vml_Matrix3 representation
 	///	@var	vml_Transform 
 	///	@var	rotation
 	///
 	typedef struct xvml_Transform
	{
 		vml_Vector3 translation;
 		vml_Matrix3 rotation;
 	} vml_Transform;
/// @}



//	@ingroup vml_Plane
/// @{
 	/// @brief	4x4 homogeneous vml_Transform vml_Matrix3 representation
 	///	@var	point 
 	///	@var	direction
 	///
 	typedef struct xvml_Plane
	{
 		vml_Vector3 point;
 		vml_Vector3 direction;
 	} vml_Plane;
/// @}








///	@ingroup vml_Vector3
/// @{

	/// @brief Gets vector X-component
	#define vml_X(vec) 	vec[0UL]

	/// @brief Gets vector Y-component
	#define vml_Y(vec) 	vec[1UL]

	/// @brief Gets vector Z-component
	#define vml_Z(vec) 	vec[2UL]


	VML_API void 	vml_VectorSetEq( vmlTy_p vecOUT, vmlTy_p vecIN, VSS_1 );
	VML_API void 	vml_VectorSetEqD2S( vmlTy_p vecOUT, vmlTy_dp vecIN, VSS_1 );
	VML_API void 	vml_VectorNegate( vmlTy_p vecOUT, vmlTy_p vecIN, VSS_1 );
	VML_API void 	vml_VectorAdd( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API void 	vml_VectorAddIPL( vmlTy_p vecOUT, vmlTy_p vecINA, VSS_1 );
	VML_API void 	vml_VectorScaleAdd( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB, rF32* scale, VSS_1 );
	VML_API void 	vml_VectorScaleAddIPL( vmlTy_p vecOUT, vmlTy_p vecINA, rF32* scale, VSS_1 );
	VML_API void 	vml_VectorSub( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API void	vml_VectorSubIPL( vmlTy_p vecOUT, vmlTy_p vecINA, VSS_1 );
	VML_API void 	vml_VectorScaleSub( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB, rF32* scale, VSS_1 );
	VML_API void 	vml_VectorScaleSubIPL( vmlTy_p vecOUT, vmlTy_p vecINA, rF32* scale, VSS_1 );
	VML_API void 	vml_VectorEMult( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API void 	vml_VectorEDiv( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API void 	vml_VectorESig( vmlTy_p vecOUT, vmlTy_p vecINA, VSS_1 );
	VML_API void 	vml_VectorScale( vmlTy_p vecOUT, vmlTy_p vecIN, rF32 *scale, VSS_1 ); 
	VML_API void 	vml_VectorScaleIPL( vmlTy_p vecOUT, rF32 *scale, VSS_1 );
	VML_API void 	vml_VectorScaleNR( vmlTy_p vecOUT, vmlTy_p vecIN, rF32 scale, VSS_1 );
	VML_API void 	vml_VectorScaleCNR( vmlTy_p vecOUT, rF32 scale, VSS_1 );
	VML_API void 	vml_VectorScaleDiv( vmlTy_p vecOUT, vmlTy_p vecIN, rF32 *scale, VSS_1 ); 
	VML_API void 	vml_VectorShift( vmlTy_p vecOUT, vmlTy_p vecIN, rF32 shift, VSS_1 );
	VML_API void 	vml_VectorShiftIPL( vmlTy_p vecOUT, rF32 shift, VSS_1 );	
	VML_API rF32 	vml_VectorDot( vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API rF32 	vml_VectorAllignment( vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API rF32 	vml_VectorNorm( vmlTy_p vecIN,  rF32 order, VSS_1 );
	VML_API void 	vml_VectorProjection( vmlTy_p vecOUT, vmlTy_p vecIN, vml_Plane* plane, VSS_1 );
	VML_API rF32 	vml_VectorL1Diff( vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API rF32	vml_VectorL2Diff( vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API rF32 	vml_VectorNormalize( vmlTy_p vecOUT, vmlTy_p vecIN, VSS_1 );
	VML_API rF32 	vml_VectorNormalizeIPL( vmlTy_p vecOUT, VSS_1 );
	VML_API rF32 	vml_VectorNormalizeDiff( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB, VSS_1 );
	VML_API void 	vml_VectorMaskMax( vmlTy_p vecOUT, vmlTy_p vecIN, VSS_1 );
	VML_API void 	vml_VectorMaskMaxIPL( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_VectorMaskMin( vmlTy_p vecOUT, vmlTy_p vecIN, VSS_1 );
	VML_API void 	vml_VectorMaskMinIPL( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_VectorCross( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB , VSS_1 );
	VML_API rF32 	vml_VectorTriple( vmlTy_p vecINA, vmlTy_p vecINB, vmlTy_p vecINC, VSS_1 );
	VML_API void 	vml_VectorAverage( vmlTy_p vecOUT, vmlTy_p vecINA, vmlTy_p vecINB , VSS_1 );
	VML_API void 	vml_VectorAbs( vmlTy_p vecOUT, vmlTy_p vecINA , VSS_1 );
	VML_API void 	vml_VectorAbsIPL( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_VectorXRotation( vmlTy_p vecOUT, vmlTy_p vecIN, rF32 *angle, VSS_1 );
	VML_API void 	vml_VectorYRotation( vmlTy_p vecOUT, vmlTy_p vecIN, rF32 *angle, VSS_1 );
	VML_API void 	vml_VectorZRotation( vmlTy_p vecOUT, vmlTy_p vecIN, rF32 *angle, VSS_1 );
	VML_API void 	vml_VectorXRotationIPL( vmlTy_p vecOUT, rF32 *angle, VSS_1 );
	VML_API void 	vml_VectorYRotationIPL( vmlTy_p vecOUT, rF32 *angle, VSS_1 );
	VML_API void 	vml_VectorZRotationIPL( vmlTy_p vecOUT, rF32 *angle, VSS_1 );
	VML_API void 	vml_VectorLoadZeros( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_VectorLoadPolar( vmlTy_p vecOUT, rF32* theta, rF32* alpha, VSS_1 );
	VML_API void 	vml_VectorLoadPolar90a( vmlTy_p vecOUT, rF32* theta, rF32* alpha, VSS_1 );
	VML_API void 	vml_VectorLoadPolar90t( vmlTy_p vecOUT, rF32* theta, rF32* alpha, VSS_1 );
	VML_API void 	vml_VectorLoadXAxis( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_VectorLoadYAxis( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_VectorLoadZAxis( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_VectorLoadAll( vmlTy_p vecOUT, rF32* num, VSS_1 );
	VML_API void 	vml_VectorLoadRand( vmlTy_p vecOUT, VSS_1 );
	VML_API void 	vml_Vector3RPY( vml_Vector3 vecOUT, vml_Matrix3 rotation );

/// @}










/// @ingroup vml_Matrix3
/// @{


	/// @brief vml_Matrix indexing macros
	#define MatXIDX(vml_Matrix,nrow,row,col) 			vml_Matrix[nrow*row + col]
	#define Mat2IDX(vml_Matrix,row,col) 				vml_Matrix[2UL*row + col]
	#define Mat3IDX(vml_Matrix,row,col) 				vml_Matrix[3UL*row + col]
	#define Mat4IDX(vml_Matrix,row,col) 				vml_Matrix[4UL*row + col]
	#define Mat5IDX(vml_Matrix,row,col) 				vml_Matrix[5UL*row + col]
	#define Mat6IDX(vml_Matrix,row,col) 				vml_Matrix[6UL*row + col]

	#define IDX_s(idx,n)								((idx<n)?idx:0UL)
	#define MatXIDX_s(vml_Matrix,nrow,ncol,row,col) 	vml_Matrix[nrow*IDX_s(row,nrow)+IDX_s(col,ncol)]
	#define Mat2IDX_s(vml_Matrix,row,col) 				vml_Matrix[2UL*IDX_s(row,2UL)+IDX_s(col,2UL)]
	#define Mat3IDX_s(vml_Matrix,row,col) 				vml_Matrix[3UL*IDX_s(row,3UL)+IDX_s(col,3UL)]
	#define Mat4IDX_s(vml_Matrix,row,col) 				vml_Matrix[4UL*IDX_s(row,4UL)+IDX_s(col,4UL)]
	#define Mat5IDX_s(vml_Matrix,row,col) 				vml_Matrix[5UL*IDX_s(row,5UL)+IDX_s(col,5UL)]
	#define Mat6IDX_s(vml_Matrix,row,col) 				vml_Matrix[6UL*IDX_s(row,6UL)+IDX_s(col,6UL)]

	/// @brief vml_Matrix3 element macros
	#define Mat3_11(vml_Matrix) 						Mat3IDX(vml_Matrix,0,0) 
	#define Mat3_12(vml_Matrix) 						Mat3IDX(vml_Matrix,0,1) 
	#define Mat3_13(vml_Matrix) 						Mat3IDX(vml_Matrix,0,2) 
	#define Mat3_21(vml_Matrix) 						Mat3IDX(vml_Matrix,1,0) 
	#define Mat3_22(vml_Matrix) 						Mat3IDX(vml_Matrix,1,1) 
	#define Mat3_23(vml_Matrix) 						Mat3IDX(vml_Matrix,1,2) 
	#define Mat3_31(vml_Matrix) 						Mat3IDX(vml_Matrix,2,0) 
	#define Mat3_32(vml_Matrix) 						Mat3IDX(vml_Matrix,2,1) 
	#define Mat3_33(vml_Matrix) 						Mat3IDX(vml_Matrix,2,2) 


	/// @brief	O = M
	VML_API void 	vml_MatrixSetEq( vmlTy_p matOUT, vmlTy_p matIN, MSS_1 );
	VML_API void 	vml_MatrixScale( vmlTy_p matOUT, vmlTy_p matIN, rF32 *scale, MSS_1 );
	VML_API void 	vml_MatrixTranspose( vmlTy_p matOUT, vmlTy_p matIN, MSS_1 );
	VML_API void 	vml_MV_Mult( vmlTy_p vecOUT, vmlTy_p matIN, vmlTy_p vecIN, MSS_1 );
	VML_API void 	vml_VM_Mult( vmlTy_p vecOUT, vmlTy_p vecIN, vmlTy_p matIN, MSS_1 );
	VML_API void 	vml_MM_Mult(  vmlTy_p matOUT, vmlTy_p matINA, vmlTy_p matINB, MSS_1, MSS_2 );
	VML_API void 	vml_MatrixLoadIdentity( vmlTy_p matOUT, MSS_1 );
	VML_API void 	vml_MatrixLoadZeros( vmlTy_p matOUT, MSS_1 );
	VML_API void 	vml_MatrixLoadOnes( vmlTy_p matOUT, MSS_1 );
	VML_API void	vml_MatrixLoadRand( vmlTy_p matOUT, MSS_1 );
	VML_API void 	vml_MatrixAdd( vmlTy_p matOUT, vmlTy_p matINA, vmlTy_p matINB, MSS_1 );
	VML_API void 	vml_MatrixSub( vmlTy_p matOUT, vmlTy_p matINA, vmlTy_p matINB, MSS_1 );
	VML_API void 	vml_MatrixAddIPL( vmlTy_p matOUT, vmlTy_p matIN, MSS_1 );
	VML_API void 	vml_MatrixSubIPL( vmlTy_p matOUT, vmlTy_p matIN, MSS_1 );
	VML_API void 	vml_MatrixScaleAddIPL( vmlTy_p matOUT, vmlTy_p matIN, rF32* scale, MSS_1 );
	VML_API void 	vml_MatrixScaleSubIPL( vmlTy_p matOUT, vmlTy_p matIN, rF32* scale, MSS_1 );
	
	VML_API void 	vml_Matrix3LoadSkew( vml_Matrix3 matOUT, vml_Vector3 vecIN );
	VML_API void 	vml_Matrix3LoadYRot( vml_Matrix3 matOUT, rF32 angle );
	VML_API void 	vml_Matrix3LoadXRot( vml_Matrix3 matOUT, rF32 angle );
	VML_API void 	vml_Matrix3LoadZRot( vml_Matrix3 matOUT, rF32 angle );
	VML_API void 	vml_Matrix3LoadRPY( vml_Matrix3 matOUT, vml_Vector3 vecRPY );

/// @}



///	@ingroup rIP
/// @{
	VML_API rBool 	vml_IPInUse();
	VML_API vmlTy_p	vml_IPGet();
	VML_API vmlTy_p	vml_IPClear();
/// @}



///	@ingroup rIPvml_Matrix
/// @brief 	 In-place operations (uses dynamic memory; not thread-safe)	
/// @{

	//VML_API vmlTy_p 	_rIPvml_MatrixTranspose( vmlTy_p matIN, MSS_1 );

/// @}



///	@ingroup vml_Transform
/// @{

	VML_API void vml_TransformSetEq( vml_Transform *homOUT,  vml_Transform *homIN );
	VML_API void vml_TransformLoadIdentity( vml_Transform *homOUT );
	VML_API void vml_TransformLoadDHParam( vml_Transform *homOUT, rF32 xoff, rF32 alpha, rF32 theta, rF32 zoff );
	VML_API void rHV_Mult( vml_Vector3 vecOUT, vml_Transform *homIN, vml_Vector3 vecIN );
	VML_API void rIHV_Mult( vml_Vector3 vecOUT, vml_Transform *homIN, vml_Vector3 vecIN  );
	VML_API void rHH_Mult( vml_Transform *homOUT, vml_Transform *homINA, vml_Transform *homINB );
	VML_API void rIHH_Mult( vml_Transform *homOUT, vml_Transform *homINA, vml_Transform *homINB );	
	VML_API void rIHH_Mult_NIPL( vml_Transform *homOUT, vml_Transform *homINA, vml_Transform *homINB );	

/// @}
	


#ifdef 	VML_WITH_OPENGL
VML_API void vml_Vector3GLRenderAsVertex( vml_Vector3 vec );
VML_API void vml_Vector3GLRenderAsNormal( vml_Vector3 vec );
#endif




///	@ingroup Defined_Constants
/// @{
#ifdef 	VML_WITH_CONSTANTS
extern 	vml_Vector3 		vml_Zeros3;
extern 	vml_Vector3 		vml_Ones3;
extern 	vml_Vector3 		vml_Xhat3;
extern 	vml_Vector3 		vml_Yhat3;
extern 	vml_Vector3 		vml_Zhat3;
extern 	vml_Vector3 		vml_Ones3n;
extern 	vml_Vector3 		vml_Xhat3n;
extern 	vml_Vector3 		vml_Yhat3n;
extern 	vml_Vector3 		vml_Zhat3n;

extern 	vml_Matrix3 		vml_I3x3;
extern 	vml_Matrix3 		vml_Zeros3x3;

extern  vml_Plane			vml_XY;
extern  vml_Plane			vml_XZ;
extern  vml_Plane			vml_YZ;
#endif
/// @}




/// @ingroup Other_Useful_Bits
/// @{

#define FILL2(num)	 	{ num, num }
#define FILL3(num)		{ num, num, num }
#define FILL4(num)		{ num, num, num, num }
#define FILL5(num)		{ num, num, num, num, num }
#define FILL6(num)		{ num, num, num, num, num, num }
#define FILL7(num)	 	{ num, num, num, num, num, num, num }
#define FILL8(num)	 	{ num, num, num, num, num, num, num, num }
#define FILL9(num)	 	{ num, num, num, num, num, num, num, num, num }
#define FILL10(num)	 	{ num, num, num, num, num, num, num, num, num, num }


#if __cplusplus
/// vmlTy_p Stream Operator
std::ostream& operator| ( std::ostream& os, const vml_Size* fmt );
std::ostream& operator| ( std::ostream& os, const unsigned int* fmt );
std::ostream& operator<<( std::ostream& os, vmlTy_p usep );
#endif

/// @}


#endif
