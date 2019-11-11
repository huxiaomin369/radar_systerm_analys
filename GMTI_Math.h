//////////////////////////////////////////////////
//定义基本的矩阵与复数运算、FFT等///
/////////////////////////////////////////////////
//   $pyHe, The RSP.
//   $Revision: 1.0 $  $Date: 2019/11/6 17:40 $
/////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

class	_declspec(dllexport)Complex;
class	_declspec(dllexport)Complex_Exp;
class	_declspec(dllexport)Complex_Vec;
class	_declspec(dllexport)Complex_Mat;


////////////////////////////////////

////////////////////////////////////
//定义复数、复矢量和复矩阵类型
class	_declspec(dllexport)Complex
{
public:
	float re;//实部
	float im;//虚部
	Complex();
	Complex(float re,float im);
	void Show();
};//x=re+j*im;


class	_declspec(dllexport)Complex_Exp
{
public:
	float rho;//模值
	float theta;//相角(rad).[-pi pi)
};//x=rho*exp(j*theta)

class	_declspec(dllexport)Complex_Vec
{
public:
	Complex* data;//
	int len;//复向量的有效长度

	Complex_Vec(int max_size);//构造函数
	~Complex_Vec();//析构函数
	void Show();
};//data[len]

class	_declspec(dllexport)Complex_Mat
{
public:
	Complex** data;//
	int nrow;//复矩阵的行数(range_num)
	int mcol;//复矩阵的列数(azimuth_num)

	Complex_Mat(int nrow,int mcol);//构造一般矩阵
	Complex_Mat(int nrow);//构造方阵
	~Complex_Mat();
	void Show();//
};//data[nrow][mcol]

#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)

////////////////
//定义内联函数
////////////////
extern inline Complex _declspec(dllexport)cAdd(Complex x,Complex y);//x+y
extern inline Complex _declspec(dllexport)cSub(Complex x,Complex y);//x-y
extern inline Complex _declspec(dllexport)cMult(Complex x,Complex y);//x*y
extern inline Complex _declspec(dllexport)cDiv(Complex x,Complex y);//x/y
extern inline Complex _declspec(dllexport)CompExp2Comp(Complex_Exp x);//复指数转成复数
extern inline float _declspec(dllexport)Modulu(Complex x);//返回复数的模值
extern inline Complex _declspec(dllexport)Conj(Complex x);//返回复数的共轭值
extern inline int _declspec(dllexport)ceil(int x,int y);//对x/y的结果向上取整.ceil(5,4)=2
extern inline int _declspec(dllexport)fix(int x,int y);//对x/y的结果向下取整.fix(5,2)=2



//////////////////////////
//指针数组分配和释放函数
//////////////////////////
void _declspec(dllexport)Init_pArray(Complex** pArray,int n,int m);
//分配具有n个元素的指针数组且指针数组的每个元素指向具有m个元素的一维数组
//pArray=new Complex* [n]自身应该在主调函数中申请

void _declspec(dllexport)Init_pArray(int** pArray,int n,int m);

void _declspec(dllexport)Discard_pArray(Complex** pArray,int n);
//释放n个元素的指针数组pArray的每一个元素(对应为每个一维数组)
//pArray自身的释放应该在主调函数中完成

void _declspec(dllexport)Discard_pArray(int** pArray,int n);

////////////////////////////////////////
//复数转换为复指数和取相角操作(-PI-PI)
////////////////////////////////////////
Complex_Exp _declspec(dllexport)Comp2CompExp(Complex x);//复数转成复指数
float _declspec(dllexport)Angle(Complex x);//返回复数的相角(rad.[-pi:pi))


////////////////////////////////////////////////////
//FFT/IFFT/FFT2/IFFT2/FFTSHIFT/IFFTSHIFT/Conv/Corr
///////////////////////////////////////////////////

int _declspec(dllexport)log2(int x);//向上取整.log2(5)=3; log2(4)=2.
//2的对数向上取整

int _declspec(dllexport)Bit_Inv(int x,int bits);//二进制逆序
//输出为输入的二进制逆序数值(用于FFT).按bits位求逆
//例如：Bit_Inv(6,5):6->00110->01100=12

void _declspec(dllexport)Pos_Inv(const Complex* in_ptr, Complex* out_ptr, int len);
//对长度为len的输入数组元素按下标进行二进制逆序重排得到新的序列out_ptr 
//子函数中不能修改in_ptr所指内容(用于FFT)

void _declspec(dllexport)FFT(const Complex* in_ptr, Complex* out_ptr, int len, int plen);
//快速傅立叶变换(FFT)
//对len点的输入序列进行plen点FFT运算,输出由out_ptr指出且长度为plen
//如果plen<len则对输入序列进行截取；如果plen>len则对输入序列添零延长
//plen 为2的整数次幂(plen=2^n)则采用快速算法，否则按照DFT计算

void _declspec(dllexport)FFT(const Complex* in_ptr, Complex* out_ptr, int len);
//对len点的输入序列进行len点FFT运算,输出由out_ptr指出且与输入等长
//len为2的整数次幂(len=2^n)则采用时域基2抽取快速算法，否则按照DFT计算

void _declspec(dllexport)FFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);  
//对in_ptr所指的数据作二维FFT并把结果保存到out_ptr中

void _declspec(dllexport)IFFT(const Complex* in_ptr, Complex* out_ptr,int len);
//快速逆傅立叶变换(IFFT)
//对len点的输入序列进行len点IFFT运算,输出由out_ptr指出且与输入等长
//len为2的整数次幂(len=2^n)则采用时域基2抽取快速算法，否则按照IDFT计算

void _declspec(dllexport)IFFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);     
//对in_ptr所指的数据作二维IFFT并把结果保存到out_ptr中


void _declspec(dllexport)FFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len);
//对len点输入序列in_ptr交换前后两部分后存入out_ptr
void _declspec(dllexport)FFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr);
//对向量in_ptr作前后(或上下)两部分交换存入向量out_ptr中
void _declspec(dllexport)FFTSHIFT(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);
//对矩阵in_ptr作前后、上下两部分交换存入矩阵out_ptr中
void _declspec(dllexport)FFTSHIFT(Complex_Mat* io_ptr);
//交换后的结果同时存入输入矩阵

void _declspec(dllexport)IFFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len);
//对len点输入序列in_ptr交换前后两部分后存入out_ptr
//IFFTSHIFT(FFTSHIFT(X))=X;
void _declspec(dllexport)IFFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr);
//向量IFFTSHIFT
void _declspec(dllexport)IFFTSHIFT(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);
//矩阵IFFTSHIFT
void _declspec(dllexport)IFFTSHIFT(Complex_Mat* io_ptr);
//交换后的结果同时存入输入矩阵

void _declspec(dllexport)Conv(const Complex_Vec* x,const Complex_Vec* h,Complex* y);
//计算输入序列x和h的线形卷积(基于FFT和IFFT)
//输出序列y=x*h(长度为length(x)+length(h)-1)
//补零后利用FFT和IFFT计算线形卷积

double _declspec(dllexport)Corr(const Complex_Vec* x,const Complex_Vec* h,Complex* y);
//计算输入序列x和h的循环相关(基于FFT和IFFT)
//返回h相对于x的偏移量(最大值对应的位置)和循环相关序列y
//输出序列y(m)=x(n).conj(h(n-m))--h(n)向右循环移位
//则y(m)=conj(IFFT(FFT(x,len).*conj(FFT(h,len))))
//输出序列长度len为x和h的最大长度值.m=[0:len-1]
//计算过程中对较短序列补零至较长序列的长度


/////////////////
//矩阵操作函数
////////////////
void _declspec(dllexport)Matrix_Mult(const Complex_Mat* pA,const Complex_Mat* pB,Complex_Mat* pC);
//矩阵pC=pA*pB

void _declspec(dllexport)Matrix_cTrans(const Complex_Mat* pA,Complex_Mat* pB);
//矩阵的共轭转置(Conjugate Transpose).pB=pA'


void _declspec(dllexport)Matrix_Trans(const Complex_Mat* pA,Complex_Mat* pB);
//矩阵的转置(Transpose).pB=pA.'

void _declspec(dllexport)Mat2Vec(const Complex_Mat* pA,Complex_Mat* pB);
//将矩阵pA以列优先转换为列矢量(即pB->mcol=1)

void _declspec(dllexport)Matrix_Inv(Complex** pRx,Complex** inv_Rx,int nrow);
//对方阵pRx[nrow][nrow]求逆并存入inv_Rx
void _declspec(dllexport)Matrix_Inv(const Complex_Mat* pA,Complex_Mat* inv_A);
//对方阵pA[nrow][nrow]求逆并存入inv_A

void _declspec(dllexport)Matrix_pInv(Complex** pRx,Complex** pinv_Rx,int nrow,int mcol);
//求矩阵pRx[nrow][mcol]的伪逆并存入pinv_Rx
//如果nrow>mcol=>pinv_Rx=inv(pRx'*pRx)*pRx'==>pinv_Rx[mcol][nrow]==>pinv_Rx*pRx=I;
//如果nrow<mcol=>pinv_Rx=pRx'*inv(pRx*pRx')==>pinv_Rx[mcol][nrow]==>pRx*pinv_Rx=I;
//满足:dim(I)=min(nrow,mcol)
void _declspec(dllexport)Matrix_pInv(const Complex_Mat* pA,Complex_Mat* pinv_A);


// 统计函数

double _declspec(dllexport)mean(const sqList* sq);
//返回数据的均值
double _declspec(dllexport)mean(const double* pd,int n);
float _declspec(dllexport)mean(const float* pd,int n);
Complex _declspec(dllexport)mean(const Complex* pd,int n);
//返回n点数组pd的均值


void _declspec(dllexport)sort(const double* pd,double* ps,int* index,int n);
//对n点数组pd的元素从小到大升序排序
//返回排序后的数组ps和对应元素在pd中的下标

float _declspec(dllexport)Find_Max_Value(const Complex_Mat* in_ptr,int* index);
//返回矩阵元素模值的最大值和对应下标
//index[0]存放行下标,index[1]存放列下标
float _declspec(dllexport)Find_Max_Value(const Complex_Vec* in_ptr,int& index);
//返回向量元素模值的最大值和对应下标(存放在index中)

void _declspec(dllexport)Ls_Fit(const Complex_Mat* A, const Complex_Mat* b, Complex_Vec* coeff);  
//最小二乘拟合  
//A*coeff=b
void _declspec(dllexoprt)Ls_Fit(const double A, const double  b, double  k);  
//实数最小二乘拟合  //于GMTI_Math中定义
//A*k=b