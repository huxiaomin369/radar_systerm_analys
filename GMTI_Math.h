//////////////////////////////////////////////////
//��������ľ����븴�����㡢FFT��///
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
//���帴������ʸ���͸���������
class	_declspec(dllexport)Complex
{
public:
	float re;//ʵ��
	float im;//�鲿
	Complex();
	Complex(float re,float im);
	void Show();
};//x=re+j*im;


class	_declspec(dllexport)Complex_Exp
{
public:
	float rho;//ģֵ
	float theta;//���(rad).[-pi pi)
};//x=rho*exp(j*theta)

class	_declspec(dllexport)Complex_Vec
{
public:
	Complex* data;//
	int len;//����������Ч����

	Complex_Vec(int max_size);//���캯��
	~Complex_Vec();//��������
	void Show();
};//data[len]

class	_declspec(dllexport)Complex_Mat
{
public:
	Complex** data;//
	int nrow;//�����������(range_num)
	int mcol;//�����������(azimuth_num)

	Complex_Mat(int nrow,int mcol);//����һ�����
	Complex_Mat(int nrow);//���췽��
	~Complex_Mat();
	void Show();//
};//data[nrow][mcol]

#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)

////////////////
//������������
////////////////
extern inline Complex _declspec(dllexport)cAdd(Complex x,Complex y);//x+y
extern inline Complex _declspec(dllexport)cSub(Complex x,Complex y);//x-y
extern inline Complex _declspec(dllexport)cMult(Complex x,Complex y);//x*y
extern inline Complex _declspec(dllexport)cDiv(Complex x,Complex y);//x/y
extern inline Complex _declspec(dllexport)CompExp2Comp(Complex_Exp x);//��ָ��ת�ɸ���
extern inline float _declspec(dllexport)Modulu(Complex x);//���ظ�����ģֵ
extern inline Complex _declspec(dllexport)Conj(Complex x);//���ظ����Ĺ���ֵ
extern inline int _declspec(dllexport)ceil(int x,int y);//��x/y�Ľ������ȡ��.ceil(5,4)=2
extern inline int _declspec(dllexport)fix(int x,int y);//��x/y�Ľ������ȡ��.fix(5,2)=2



//////////////////////////
//ָ�����������ͷź���
//////////////////////////
void _declspec(dllexport)Init_pArray(Complex** pArray,int n,int m);
//�������n��Ԫ�ص�ָ��������ָ�������ÿ��Ԫ��ָ�����m��Ԫ�ص�һά����
//pArray=new Complex* [n]����Ӧ������������������

void _declspec(dllexport)Init_pArray(int** pArray,int n,int m);

void _declspec(dllexport)Discard_pArray(Complex** pArray,int n);
//�ͷ�n��Ԫ�ص�ָ������pArray��ÿһ��Ԫ��(��ӦΪÿ��һά����)
//pArray������ͷ�Ӧ�����������������

void _declspec(dllexport)Discard_pArray(int** pArray,int n);

////////////////////////////////////////
//����ת��Ϊ��ָ����ȡ��ǲ���(-PI-PI)
////////////////////////////////////////
Complex_Exp _declspec(dllexport)Comp2CompExp(Complex x);//����ת�ɸ�ָ��
float _declspec(dllexport)Angle(Complex x);//���ظ��������(rad.[-pi:pi))


////////////////////////////////////////////////////
//FFT/IFFT/FFT2/IFFT2/FFTSHIFT/IFFTSHIFT/Conv/Corr
///////////////////////////////////////////////////

int _declspec(dllexport)log2(int x);//����ȡ��.log2(5)=3; log2(4)=2.
//2�Ķ�������ȡ��

int _declspec(dllexport)Bit_Inv(int x,int bits);//����������
//���Ϊ����Ķ�����������ֵ(����FFT).��bitsλ����
//���磺Bit_Inv(6,5):6->00110->01100=12

void _declspec(dllexport)Pos_Inv(const Complex* in_ptr, Complex* out_ptr, int len);
//�Գ���Ϊlen����������Ԫ�ذ��±���ж������������ŵõ��µ�����out_ptr 
//�Ӻ����в����޸�in_ptr��ָ����(����FFT)

void _declspec(dllexport)FFT(const Complex* in_ptr, Complex* out_ptr, int len, int plen);
//���ٸ���Ҷ�任(FFT)
//��len����������н���plen��FFT����,�����out_ptrָ���ҳ���Ϊplen
//���plen<len����������н��н�ȡ�����plen>len����������������ӳ�
//plen Ϊ2����������(plen=2^n)����ÿ����㷨��������DFT����

void _declspec(dllexport)FFT(const Complex* in_ptr, Complex* out_ptr, int len);
//��len����������н���len��FFT����,�����out_ptrָ����������ȳ�
//lenΪ2����������(len=2^n)�����ʱ���2��ȡ�����㷨��������DFT����

void _declspec(dllexport)FFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);  
//��in_ptr��ָ����������άFFT���ѽ�����浽out_ptr��

void _declspec(dllexport)IFFT(const Complex* in_ptr, Complex* out_ptr,int len);
//�����渵��Ҷ�任(IFFT)
//��len����������н���len��IFFT����,�����out_ptrָ����������ȳ�
//lenΪ2����������(len=2^n)�����ʱ���2��ȡ�����㷨��������IDFT����

void _declspec(dllexport)IFFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);     
//��in_ptr��ָ����������άIFFT���ѽ�����浽out_ptr��


void _declspec(dllexport)FFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len);
//��len����������in_ptr����ǰ�������ֺ����out_ptr
void _declspec(dllexport)FFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr);
//������in_ptr��ǰ��(������)�����ֽ�����������out_ptr��
void _declspec(dllexport)FFTSHIFT(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);
//�Ծ���in_ptr��ǰ�����������ֽ����������out_ptr��
void _declspec(dllexport)FFTSHIFT(Complex_Mat* io_ptr);
//������Ľ��ͬʱ�����������

void _declspec(dllexport)IFFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len);
//��len����������in_ptr����ǰ�������ֺ����out_ptr
//IFFTSHIFT(FFTSHIFT(X))=X;
void _declspec(dllexport)IFFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr);
//����IFFTSHIFT
void _declspec(dllexport)IFFTSHIFT(const Complex_Mat* in_ptr, Complex_Mat* out_ptr);
//����IFFTSHIFT
void _declspec(dllexport)IFFTSHIFT(Complex_Mat* io_ptr);
//������Ľ��ͬʱ�����������

void _declspec(dllexport)Conv(const Complex_Vec* x,const Complex_Vec* h,Complex* y);
//������������x��h�����ξ��(����FFT��IFFT)
//�������y=x*h(����Ϊlength(x)+length(h)-1)
//���������FFT��IFFT�������ξ��

double _declspec(dllexport)Corr(const Complex_Vec* x,const Complex_Vec* h,Complex* y);
//������������x��h��ѭ�����(����FFT��IFFT)
//����h�����x��ƫ����(���ֵ��Ӧ��λ��)��ѭ���������y
//�������y(m)=x(n).conj(h(n-m))--h(n)����ѭ����λ
//��y(m)=conj(IFFT(FFT(x,len).*conj(FFT(h,len))))
//������г���lenΪx��h����󳤶�ֵ.m=[0:len-1]
//��������жԽ϶����в������ϳ����еĳ���


/////////////////
//�����������
////////////////
void _declspec(dllexport)Matrix_Mult(const Complex_Mat* pA,const Complex_Mat* pB,Complex_Mat* pC);
//����pC=pA*pB

void _declspec(dllexport)Matrix_cTrans(const Complex_Mat* pA,Complex_Mat* pB);
//����Ĺ���ת��(Conjugate Transpose).pB=pA'


void _declspec(dllexport)Matrix_Trans(const Complex_Mat* pA,Complex_Mat* pB);
//�����ת��(Transpose).pB=pA.'

void _declspec(dllexport)Mat2Vec(const Complex_Mat* pA,Complex_Mat* pB);
//������pA��������ת��Ϊ��ʸ��(��pB->mcol=1)

void _declspec(dllexport)Matrix_Inv(Complex** pRx,Complex** inv_Rx,int nrow);
//�Է���pRx[nrow][nrow]���沢����inv_Rx
void _declspec(dllexport)Matrix_Inv(const Complex_Mat* pA,Complex_Mat* inv_A);
//�Է���pA[nrow][nrow]���沢����inv_A

void _declspec(dllexport)Matrix_pInv(Complex** pRx,Complex** pinv_Rx,int nrow,int mcol);
//�����pRx[nrow][mcol]��α�沢����pinv_Rx
//���nrow>mcol=>pinv_Rx=inv(pRx'*pRx)*pRx'==>pinv_Rx[mcol][nrow]==>pinv_Rx*pRx=I;
//���nrow<mcol=>pinv_Rx=pRx'*inv(pRx*pRx')==>pinv_Rx[mcol][nrow]==>pRx*pinv_Rx=I;
//����:dim(I)=min(nrow,mcol)
void _declspec(dllexport)Matrix_pInv(const Complex_Mat* pA,Complex_Mat* pinv_A);


// ͳ�ƺ���

double _declspec(dllexport)mean(const sqList* sq);
//�������ݵľ�ֵ
double _declspec(dllexport)mean(const double* pd,int n);
float _declspec(dllexport)mean(const float* pd,int n);
Complex _declspec(dllexport)mean(const Complex* pd,int n);
//����n������pd�ľ�ֵ


void _declspec(dllexport)sort(const double* pd,double* ps,int* index,int n);
//��n������pd��Ԫ�ش�С������������
//��������������ps�Ͷ�ӦԪ����pd�е��±�

float _declspec(dllexport)Find_Max_Value(const Complex_Mat* in_ptr,int* index);
//���ؾ���Ԫ��ģֵ�����ֵ�Ͷ�Ӧ�±�
//index[0]������±�,index[1]������±�
float _declspec(dllexport)Find_Max_Value(const Complex_Vec* in_ptr,int& index);
//��������Ԫ��ģֵ�����ֵ�Ͷ�Ӧ�±�(�����index��)

void _declspec(dllexport)Ls_Fit(const Complex_Mat* A, const Complex_Mat* b, Complex_Vec* coeff);  
//��С�������  
//A*coeff=b
void _declspec(dllexoprt)Ls_Fit(const double A, const double  b, double  k);  
//ʵ����С�������  //��GMTI_Math�ж���
//A*k=b