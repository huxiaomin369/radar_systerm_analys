#include "CFAR_SubFunction.h"


/////////////////////////////////////////////////
//实现GMTI数据仿真和内部处理数据格式
/////////////////////////////////////////////////
//   $zwYang, The RSP.
//   $Revision: 1.0 $  $Date: 2009/04/14 18:20 $
/////////////////////////////////////////////////
#define eb_eps 1e-5 //基线估计的内部精度
#define TRUE 1
#define FALSE 0

////////////////////////////
//GMTI数据处理内部数据类型
///////////////////////////
//复数
Complex::Complex() { re=0; im=0;}
Complex::Complex(float re,float im)
{ this->re=re; this->im=im;}
void Complex::Show()
{ cout<<re<<"+i"<<im<<endl;}


//复矢量
Complex_Vec::Complex_Vec(int max_size)
{
	if(max_size > 0)
	{
		data=new Complex[max_size];	assert(data!=NULL);//条件不成立则终止
		len=max_size;
		for(int i=0;i<max_size;i++) { data[i].re=0; data[i].im=0;}
	}
}//
Complex_Vec::~Complex_Vec()
{
	delete[] data;
	data=NULL;
}//
void Complex_Vec::Show()
{
	for (int i=0;i<len;i++){ cout<<data[i].re<<"+i"<<data[i].im<<" ";}
	cout<<endl;
}//

//复矩阵
Complex_Mat::Complex_Mat(int nrow,int mcol)	
{
	assert((nrow > 0) && (mcol > 0));
	if((nrow > 0) && (mcol > 0))
	{
		data=new Complex* [nrow];	assert(data!=NULL);//point array 
		for(int n=0; n<nrow; n++)
		{	
			data[n]=new Complex[mcol]; assert(data[n]!=NULL);//
			for(int m=0;m<mcol;m++) {	data[n][m].re=0; data[n][m].im=0;}
		}
		this->nrow=nrow;
		this->mcol=mcol;
	}
}//构造一般矩阵.data[nrow][mcol]
Complex_Mat::Complex_Mat(int nrow)	
{
	assert(nrow > 0);
	if(nrow > 0)
	{
		data=new Complex* [nrow]; assert(data!=NULL);//point array 
		for(int n=0; n<nrow; n++)
		{	
			data[n]=new Complex[nrow]; assert(data[n]!=NULL);//
			for(int m=0;m<nrow;m++) {	data[n][m].re=0; data[n][m].im=0;}
		}
		this->nrow=nrow;
		this->mcol=nrow;
	}
}//构造方阵

Complex_Mat::~Complex_Mat()
{
	for(int n=0;n<nrow;n++)
	{
		delete[] data[n];
		data[n]=NULL;
	}
	delete[] data; data=NULL;
}//析构函数

void Complex_Mat::Show()
{
	for (int i=0;i<nrow;i++)
	{
		for (int j=0;j<mcol;j++)
		{
			cout<<data[i][j].re<<"+i"<<data[i][j].im<<" ";
		}
		cout<<endl;
	}
}

/////////////////////////////////////////////////
//实现公共的数字信号处理子函数
/////////////////////////////////////////////////
//   $zwYang, The RSP.
//   $Revision: 1.0 $  $Date: 2009/04/16 17:20 $
/////////////////////////////////////////////////

////////////////
//实现内联函数
////////////////
extern inline Complex cAdd(Complex x,Complex y)
{
	Complex r;
	r.re=x.re+y.re;
	r.im=x.im+y.im;
	return r;
}//x+y

extern inline Complex cSub(Complex x,Complex y)
{
	Complex r;
	r.re=x.re-y.re;
	r.im=x.im-y.im;
	return r;
}//x-y

extern inline Complex cMult(Complex x,Complex y)
{
	Complex r;
	r.re=x.re*y.re-x.im*y.im;
	r.im=x.re*y.im+x.im*y.re;
	return r;
}//x*y

extern inline Complex cDiv(Complex x,Complex y)
{
    Complex r;
	double ay;//ay=(Modulu(y)).^2 
	ay=(y.re*y.re+y.im*y.im);
	r.re=(float) ((x.re*y.re+x.im*y.im)/ay);
	r.im=(float) ((x.im*y.re-x.re*y.im)/ay);
	return r;
}//x/y

//复指数转成复数
extern inline Complex CompExp2Comp(Complex_Exp x)
{
	Complex y;
	y.re=(float) (x.rho*cos(x.theta));
	y.im=(float) (x.rho*sin(x.theta));
	return y;
}

//返回复数的模值
extern inline float Modulu(Complex x)
{
	float r=(float) (sqrt(x.re*x.re+x.im*x.im));
	return r;
}

//返回复数的共轭值
extern inline Complex Conj(Complex x)
{
	Complex y;
	y.re=x.re; y.im=(-1)*x.im;
	return y;
}

//对x/y的结果向上取整.ceil(5,4)=2
extern inline int ceil(int x,int y)
{
	int s=x/y;
	if((s*y) < x){	s=s+1;}
	return s;
}

//对x/y的结果向下取整.fix(5,2)=2
extern inline int fix(int x,int y)
{
	int s=x/y;
	if((s*y) > x){	s=s-1;}
	return s;
}

//////////////////////
//指针数组分配与释放 
//////////////////////
//分配具有n个元素的指针数组且指针数组的每个元素指向具有m个元素的一维数组
//pArray=new Complex* [n]自身应该在主调函数中申请
void Init_pArray(Complex** pArray,int n,int m)
{
	for(int i=0;i<n;i++)
	{
		pArray[i]=new Complex[m];
//		if(pArray[i]==NULL) {	exit(1);}
		assert(pArray[i]!=NULL);
	}
}

void Init_pArray(int** pArray,int n,int m)
{
	for(int i=0;i<n;i++)
	{
		pArray[i]=new int[m];
		assert(pArray[i]!=NULL);
	}
}

//释放n个元素的指针数组pArray的每一个元素(对应为每个一维数组)
//pArray自身的释放应该在主调函数中完成
void Discard_pArray(Complex** pArray,int n)
{
	for(int i=0;i<n;i++)
	{
		delete[] pArray[i];
		pArray[i]=NULL;
	}
}

void Discard_pArray(int** pArray,int n)
{
	for(int i=0;i<n;i++)
	{
		delete[] pArray[i];
		pArray[i]=NULL;
	}
}


////////////////////////////////////////
//复数转换为复指数和取相角操作(-PI-PI)
////////////////////////////////////////
//复数转成复指数
Complex_Exp Comp2CompExp(Complex x)
{
	Complex_Exp y;
	y.rho=(float) (sqrt(x.re*x.re+x.im*x.im));
	if(x.re >=0 )
	{
		if((x.re<EPS) & ((x.im<EPS) & (x.im> (-1)*EPS)))
		{	y.theta=0;}
		else
		{	y.theta=(float) atan(x.im/x.re);}
	}//[-PI/2:PI/2]
	else
	{
		if((x.re>(-1)*EPS) & ((x.im<EPS) & (x.im> (-1)*EPS)))
		{	y.theta=0;}
		else 
		{
			if(x.im < (-1)*EPS) 
			{	y.theta=(float) (atan(x.im/x.re)-PI);}//[-PI:-PI/2]
			else
			{	y.theta=(float) (atan(x.im/x.re)+PI);}//[PI/2:PI]
		}
	}
	return y;
}

//返回复数的相角(rad.[-pi:pi))
float Angle(Complex x)
{
	
	double theta; 
	int sign=1;//x.im的符号.1:表示正;-1表示负
	if(x.im<0) { sign=-1;} 
	if(x.re >=0 )
	{
		if((x.re<EPS) & ((x.im<EPS) & (x.im> (-1)*EPS)))
		{	theta=0;}
		else
		{	theta=atan(x.im/x.re);}
	}//[-PI/2:PI/2]
	else
	{
		if((x.re>(-1)*EPS) & ((x.im<EPS) & (x.im> (-1)*EPS)))
		{	theta=0;}
		else 
		{	theta=atan(x.im/x.re)+sign*PI;}
	}//[-PI:-PI/2] or [PI/2:PI]
	return (float) theta;
}

//////////////////////
//FFT/IFFT/Conv/Corr
//////////////////////
//2的对数向上取整
int log2(int x)
{
	int s=1;
	int r=0;//
	while(s<x)
	{
		s=s*2;
		r++;
	}
	return (r);
}//log2(4)=2;log2(5)=3

//输出为输入的二进制逆序数值。按bits位求逆
//例如：Bit_Inv(6,5):6->00110->01100=12
int Bit_Inv(int x,int bits)
{
	int s=0; //s为要返回的数 
	int index=bits-1; 
	while(index>=0) 
	{   
		s+=(x%2)*(int)pow(2,index); 
		x=x/2; 
		index--;
	} 
	return s; 
}//

//对输入数组元素按下标进行二进制逆序重排得到新的序列 
void Pos_Inv(const Complex* in_ptr, Complex* out_ptr, int len)
{
	int j=0;
	int bits=log2(len);
	for(int i=0;i<len;i++)
	{
		j=Bit_Inv(i,bits);
		out_ptr[j]=in_ptr[i];
	}
}


//快速傅立叶变换(FFT)
//对len点的输入序列进行plen点FFT运算,输出由out_ptr指出且长度为plen
//如果plen<len则对输入序列进行截取；如果plen>len则对输入序列添零延长
//plen 为2的整数次幂(plen=2^n)则采用快速算法，否则按照DFT计算
void FFT(const Complex* in_ptr, Complex* out_ptr, int len, int plen)
{
	Complex* inn_iptr=new Complex[plen];  assert(inn_iptr!=NULL);
	if(len >= plen)
	{
		for(int i=0;i<plen;i++)	{	inn_iptr[i]=in_ptr[i];}
	}
	else
	{
		for(int i=0;i<len;i++)	{	inn_iptr[i]=in_ptr[i];}
		for(i=len;i<plen;i++)	{	inn_iptr[i].re=0;	inn_iptr[i].im=0;}
	}//根据变换点数plen对原序列进行截取或补零

	FFT(inn_iptr,out_ptr,plen);//调用plen点FFT
    delete[] inn_iptr; inn_iptr=NULL;
}

//快速傅立叶变换
//对len点的输入序列进行len点FFT运算,输出由out_ptr指出且与输入等长
//len为2的整数次幂(len=2^n)则采用时域基2抽取快速算法，否则按照DFT计算
void FFT(const Complex* in_ptr, Complex* out_ptr, int len)
{
	int stage=log2(len);
	int flag=(int)pow(2,stage);
	if( flag==len){	flag=1;}//len=2^n
	else {flag=0;}

	double pha=((-2)*PI/len);//DFT角频率间隔(rad)
	Complex ss,dd,tt;//

	if(flag==1)
	{
		int pn=len/2;//每一级的蝶形总单元数
		Complex* W_f=new Complex[pn]; assert(W_f!=NULL);//旋转因子
		for(int i=0;i<pn;i++)
		{
			W_f[i].re=(float)cos(pha*i);
			W_f[i].im=(float)sin(pha*i);
		}//计算旋转因子

		Pos_Inv(in_ptr,out_ptr,len);

		int step=1;
		for(i=1;i<=stage;i++)
		{
			step=step*2;//step=(int)pow(2,i);
			int step_mid=step/2;
			int wskip=len/step;
			for(int k=0;k<len;k+=step)
			{
				for(int m=0;m<step_mid;m++)
				{
					int wp=m*wskip;
					int fp=k+m;//
					int bp=fp+step_mid;//

					tt=cMult(out_ptr[bp],W_f[wp]);			
					ss=cAdd(out_ptr[fp],tt);
					dd=cSub(out_ptr[fp],tt);

					out_ptr[fp]=ss;
					out_ptr[bp]=dd;
				}
			}
		}
		delete[] W_f;	W_f=NULL;
	}//len=2^n采用快速算法
	else
	{
		Complex W_f;
		for(int k=0;k<len;k++)
		{
			ss.re=0;	ss.im=0;//cummulate sum
			for(int n=0;n<len;n++)
			{
				W_f.re=(float)cos(pha*k*n);
				W_f.im=(float)sin(pha*k*n);
				tt=cMult(in_ptr[n],W_f);
				ss.re+=tt.re;  ss.im+=tt.im;
			}
			out_ptr[k]=ss;
		}
	}//DFT计算方式
}//reload FFT function


//对in_ptr作二维FFT并把结果保存到out_ptr
void FFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr)
{
	assert((in_ptr->nrow==out_ptr->nrow)&&(in_ptr->mcol==out_ptr->mcol));	
	Complex* xx=new Complex[in_ptr->nrow]; assert(xx!=NULL);//列数据
	Complex* fxx=new Complex[in_ptr->nrow]; assert(fxx!=NULL);//变换的谱

	for (int i=0;i<in_ptr->nrow;i++)
	{	FFT(in_ptr->data[i],out_ptr->data[i],in_ptr->mcol);}//沿行变换

	for (i=0;i<out_ptr->mcol;i++)
	{
		for (int j=0;j<out_ptr->nrow;j++){	xx[j]=out_ptr->data[j][i];}
		FFT(xx,fxx,in_ptr->nrow);
		for (j=0;j<out_ptr->nrow;j++){	out_ptr->data[j][i]=fxx[j];}
	}//沿列变换
	delete[] xx; xx=NULL;		delete[] fxx; fxx=NULL;
}

//快速逆傅立叶变换(IFFT)
//对len点的输入序列进行len点IFFT运算,输出由out_ptr指出且与输入等长
//len为2的整数次幂(len=2^n)则采用时域基2抽取快速算法，否则按照IDFT计算
void IFFT(const Complex* in_ptr, Complex* out_ptr,int len)
{
	int stage=log2(len);
	int flag=(int)pow(2,stage);
	if( flag==len){	flag=1;}//len=2^n
	else {flag=0;}

	double pha=((2)*PI/len);//IDFT角频率间隔(rad)
	Complex ss,dd,tt;//

	if(flag==1)
	{
		int pn=len/2;//每一级的蝶形总单元数
		Complex* W_f=new Complex[pn]; assert(W_f!=NULL);//旋转因子
		for(int i=0;i<pn;i++)
		{
			W_f[i].re=(float)cos(pha*i);
			W_f[i].im=(float)sin(pha*i);
		}//计算旋转因子

		Pos_Inv(in_ptr,out_ptr,len);

		int step=1;
		for(i=1;i<=stage;i++)
		{
			step=step*2;//step=(int)pow(2,i);
			int step_mid=step/2;
			int wskip=len/step;
			for(int k=0;k<len;k+=step)
			{
				for(int m=0;m<step_mid;m++)
				{
					int wp=m*wskip;
					int fp=k+m;//
					int bp=fp+step_mid;//

					tt=cMult(out_ptr[bp],W_f[wp]);			
					ss=cAdd(out_ptr[fp],tt);
					dd=cSub(out_ptr[fp],tt);

					out_ptr[fp].re=ss.re/2;  out_ptr[fp].im=ss.im/2;
					out_ptr[bp].re=dd.re/2;  out_ptr[bp].im=dd.im/2;
				}
			}
		}
		delete[] W_f;	W_f=NULL;
	}//len=2^n采用快速算法
	else
	{
		Complex W_f;
		Complex* inn_iptr=NULL;
		inn_iptr=new Complex[len]; assert(inn_iptr!=NULL);
		for(int n=0;n<len;n++)
		{	
			inn_iptr[n].re=in_ptr[n].re/len;
			inn_iptr[n].im=in_ptr[n].im/len;
		}

		for(n=0;n<len;n++)
		{
			ss.re=0;	ss.im=0;//cummulate sum
			for(int k=0;k<len;k++)
			{
				W_f.re=(float)cos(pha*k*n);
				W_f.im=(float)sin(pha*k*n);
				tt=cMult(inn_iptr[k],W_f);
				ss.re+=tt.re;  ss.im+=tt.im;
			}
			out_ptr[n]=ss;
		}
		delete[] inn_iptr;	inn_iptr=NULL;
	}//IDFT计算方式
}


//对in_ptr作二维IFFT并把结果保存到out_ptr
void IFFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr)
{
	assert((in_ptr->nrow==out_ptr->nrow)&&(in_ptr->mcol==out_ptr->mcol));
	Complex* xx=new Complex[in_ptr->nrow]; assert(xx!=NULL);//列数据
	Complex* fxx=new Complex[in_ptr->nrow]; assert(fxx!=NULL);//变换的谱

	for (int i=0;i<in_ptr->nrow;i++)
	{ IFFT(in_ptr->data[i],out_ptr->data[i],in_ptr->mcol);}//沿行变换

	for (i=0;i<out_ptr->mcol;i++)
	{
		for (int j=0;j<out_ptr->nrow;j++){ xx[j]=out_ptr->data[j][i];}	
		IFFT(xx,fxx,in_ptr->nrow);
		for (j=0;j<out_ptr->nrow;j++){ out_ptr->data[j][i]=fxx[j];}
	}//沿列变换
	delete[] xx; xx=NULL;		delete[] fxx; fxx=NULL;
}

//对len点输入序列in_ptr交换前后两部分后存入out_ptr
void FFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len)
{
	int pos=ceil(len,2);//前后两部分的分界点下标
	int i,j=0;
	for(i=pos;i<len;i++){ out_ptr[j]=in_ptr[i];	j++;}//后半部分
	for(i=0;i<pos;i++){ out_ptr[j]=in_ptr[i]; j++;}//前半部分
}

//对向量in_ptr作前后(或上下)两部分交换存入向量out_ptr中
void FFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr)
{
	assert(in_ptr->len==out_ptr->len);
	int pos=ceil(in_ptr->len,2);
	int i,j=0;
	for(i=pos;i<in_ptr->len;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
	for(i=0;i<pos;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
}

//对矩阵in_ptr作前后、上下两部分交换存入矩阵out_ptr中
void FFTSHIFT(const Complex_Mat* in_ptr, Complex_Mat* out_ptr)
{
	assert((in_ptr->nrow==out_ptr->nrow)&&(in_ptr->mcol==out_ptr->mcol));
	Complex* temp=new Complex[in_ptr->mcol]; assert(temp!=NULL);

	int i,j,pos;
	if(in_ptr->nrow==1) 
	{ 
		pos=ceil(in_ptr->mcol,2); j=0;
		for(i=pos;i<in_ptr->mcol;i++){ out_ptr->data[0][j]=in_ptr->data[0][i]; j++;}
		for(i=0;i<pos;i++){ out_ptr->data[0][j]=in_ptr->data[0][i]; j++;}
	}//交换列(行向量)
	else if (in_ptr->mcol==1)
	{
		pos=ceil(in_ptr->nrow,2); j=0;
		for(i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
		for(i=0;i<pos;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
	}//交换行(列向量)
	else
	{
		pos=ceil(in_ptr->nrow,2);
		for(int p=0;p<in_ptr->mcol;p++)
		{
			j=0;
			for (i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][p]=in_ptr->data[i][p]; j++;}
			for (i=0;i<pos;i++){ out_ptr->data[j][p]=in_ptr->data[i][p];j++;}	
		}//每一列前后两部分交换后的结果存入out_ptr对应列
		pos=ceil(in_ptr->mcol,2);
		for(p=0;p<in_ptr->nrow;p++)
		{
			j=0;
			for(i=pos;i<in_ptr->mcol;i++){ temp[j]=out_ptr->data[p][i]; j++;}
			for(i=0;i<pos;i++){ temp[j]=out_ptr->data[p][i]; j++;}//暂存入temp
			for(i=0;i<in_ptr->mcol;i++){ out_ptr->data[p][i]=temp[i];}//存入out_ptr
		}//每一行前后两部分交换后的结果存入out_ptr对应列
	}//
	delete[] temp; temp=NULL;
}

//对io_ptr所指的输入矩阵作前后、上下两部分交换存入输入矩阵中
void FFTSHIFT(Complex_Mat* io_ptr)
{
	int mdim=max(io_ptr->mcol,io_ptr->nrow);
	Complex* temp=new Complex[mdim]; assert(temp!=NULL);
	int i,j,pos;
	pos=ceil(io_ptr->nrow,2);
	for(int p=0;p<io_ptr->mcol;p++)
	{
		j=0;
		for (i=pos;i<io_ptr->nrow;i++){ temp[j]=io_ptr->data[i][p]; j++;}
		for (i=0;i<pos;i++){ temp[j]=io_ptr->data[i][p];j++;}
		for(i=0;i<io_ptr->nrow;i++){ io_ptr->data[i][p]=temp[i];}//存入io_ptr
	}//每一列前后两部分交换后的结果存入io_ptr对应列
	pos=ceil(io_ptr->mcol,2);
	for(p=0;p<io_ptr->nrow;p++)
	{
		j=0;
		for(i=pos;i<io_ptr->mcol;i++){ temp[j]=io_ptr->data[p][i]; j++;}
		for(i=0;i<pos;i++){ temp[j]=io_ptr->data[p][i]; j++;}//暂存入temp
		for(i=0;i<io_ptr->mcol;i++){ io_ptr->data[p][i]=temp[i];}//存入out_ptr
	}//每一行前后两部分交换后的结果存入io_ptr对应列	
	delete[] temp; temp=NULL;
}

//对len点输入序列in_ptr交换前后两部分后存入out_ptr
//IFFTSHIFT(FFTSHIFT(X))=X;
void IFFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len)
{
	int pos=fix(len,2);
	int i,j=0;
	for(i=pos;i<len;i++){ out_ptr[j]=in_ptr[i]; j++;}//后半部分
	for(i=0;i<pos;i++){ out_ptr[j]=in_ptr[i]; j++;}//前半部分
}

//向量IFFTSHIFT
void IFFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr)
{
	assert(in_ptr->len==out_ptr->len);
	int pos=fix(in_ptr->len,2);
	int i,j=0;
	for(i=pos;i<in_ptr->len;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
	for(i=0;i<pos;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
}

//矩阵IFFTSHIFT
void IFFTSHIFT(const Complex_Mat* in_ptr, Complex_Mat* out_ptr)
{
	assert((in_ptr->nrow==out_ptr->nrow)&&(in_ptr->mcol==out_ptr->mcol));
	Complex* temp=new Complex[in_ptr->mcol]; assert(temp!=NULL);

	int i,j,pos;
	if(in_ptr->nrow==1) 
	{ 
		pos=fix(in_ptr->mcol,2); j=0;
		for(i=pos;i<in_ptr->mcol;i++){ out_ptr->data[0][j]=in_ptr->data[0][i]; j++;}
		for(i=0;i<pos;i++){ out_ptr->data[0][j]=in_ptr->data[0][i]; j++;}
	}//交换列(行向量)
	else if (in_ptr->mcol==1)
	{
		pos=fix(in_ptr->nrow,2); j=0;
		for(i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
		for(i=0;i<pos;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
	}//交换行(列向量)
	else
	{
		pos=fix(in_ptr->nrow,2);
		for(int p=0;p<in_ptr->mcol;p++)
		{
			j=0;
			for (i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][p]=in_ptr->data[i][p]; j++;}
			for (i=0;i<pos;i++){ out_ptr->data[j][p]=in_ptr->data[i][p];j++;}	
		}//每一列前后两部分交换后的结果存入out_ptr对应列
		pos=fix(in_ptr->mcol,2);
		for(p=0;p<in_ptr->nrow;p++)
		{
			j=0;
			for(i=pos;i<in_ptr->mcol;i++){ temp[j]=out_ptr->data[p][i]; j++;}
			for(i=0;i<pos;i++){ temp[j]=out_ptr->data[p][i]; j++;}//暂存入temp
			for(i=0;i<in_ptr->mcol;i++){ out_ptr->data[p][i]=temp[i];}//存入out_ptr
		}//每一行前后两部分交换后的结果存入out_ptr对应列
	}//
	delete[] temp; temp=NULL;
}

//对io_ptr所指的输入矩阵作前后、上下两部分交换存入输入矩阵中
void IFFTSHIFT(Complex_Mat* io_ptr)
{
	int mdim=max(io_ptr->mcol,io_ptr->nrow);
	Complex* temp=new Complex[mdim]; assert(temp!=NULL);
	int i,j,pos;
	pos=fix(io_ptr->nrow,2);
	for(int p=0;p<io_ptr->mcol;p++)
	{
		j=0;
		for (i=pos;i<io_ptr->nrow;i++){ temp[j]=io_ptr->data[i][p]; j++;}
		for (i=0;i<pos;i++){ temp[j]=io_ptr->data[i][p];j++;}
		for(i=0;i<io_ptr->nrow;i++){ io_ptr->data[i][p]=temp[i];}//存入io_ptr
	}//每一列前后两部分交换后的结果存入io_ptr对应列
	pos=fix(io_ptr->mcol,2);
	for(p=0;p<io_ptr->nrow;p++)
	{
		j=0;
		for(i=pos;i<io_ptr->mcol;i++){ temp[j]=io_ptr->data[p][i]; j++;}
		for(i=0;i<pos;i++){ temp[j]=io_ptr->data[p][i]; j++;}//暂存入temp
		for(i=0;i<io_ptr->mcol;i++){ io_ptr->data[p][i]=temp[i];}//存入out_ptr
	}//每一行前后两部分交换后的结果存入io_ptr对应列	
	delete[] temp; temp=NULL;
}

//计算线形卷积(基于FFT和IFFT)
//输入序列x和h
//输出序列y=x*h(长度为length(x)+length(h)-1)
//补零后利用FFT和IFFT计算线形卷积
void Conv(const Complex_Vec* x,const Complex_Vec*h,Complex* y)
{
	int x_len,h_len,y_len,plen;
	x_len=x->len;	h_len=h->len;	y_len=x_len+h_len-1;
	plen=log2(y_len);	plen=(int) pow(2,plen);//FFT的变换点数

	Complex* xx=new Complex[plen];  assert(xx!=NULL);//补零后的输入序列
	Complex* hh=new Complex[plen];  assert(hh!=NULL);//补零后的输入序列
	Complex* fxx=new Complex[plen];	assert(fxx!=NULL);//FFT的变换结果
	Complex* fhh=new Complex[plen];	assert(fhh!=NULL);//FFT的变换结果
	
	//pad zeros for xx and hh
	for(int i=0;i<x_len;i++){	xx[i]=x->data[i];}
	for(i=x_len;i<plen;i++){	xx[i].re=0;	xx[i].im=0;}
	for(i=0;i<h_len;i++){	hh[i]=h->data[i];}
	for(i=h_len;i<plen;i++){	hh[i].re=0;	hh[i].im=0;}

	FFT(xx,fxx,plen);
	FFT(hh,fhh,plen);	
	delete[] xx; xx=NULL;	delete[] hh; hh=NULL;

	for(i=0;i<plen;i++){	fxx[i]=cMult(fxx[i],fhh[i]);}

	IFFT(fxx,fhh,plen);

	for(i=0;i<y_len;i++){	y[i]=fhh[i];}//截取卷积的有效值

	delete[] fxx; fxx=NULL;	 delete[] fhh; fhh=NULL;
}


//计算输入序列x和h的循环相关(基于FFT和IFFT)
//返回h相对于x的偏移量(最大值对应的位置)和循环相关序列y
//输出序列y(m)=x(n).conj(h(n-m))--h(n)向右循环移位
//      则y(m)=conj(IFFT(FFT(x,len).*conj(FFT(h,len))))
//输出序列长度len为x和h的最大长度值.m=[0:len-1]
//计算过程中对较短序列补零至较长序列的长度
double Corr(const Complex_Vec* x,const Complex_Vec* h,Complex* y)
{
	int x_len=x->len;
	int h_len=h->len;
	int plen=(x_len>h_len?x_len:h_len);

	Complex* xx=new Complex[plen];  assert(xx!=NULL);//补零后的输入序列
	Complex* hh=new Complex[plen];  assert(hh!=NULL);//补零后的输入序列
	Complex* fxx=new Complex[plen];	assert(fxx!=NULL);//FFT的变换结果
	Complex* fhh=new Complex[plen];	assert(fhh!=NULL);//FFT的变换结果
	
	//pad zeros for xx and hh
	for(int i=0;i<x_len;i++){	xx[i]=x->data[i];}
	for(i=x_len;i<plen;i++){	xx[i].re=0;	xx[i].im=0;}
	for(i=0;i<h_len;i++){	hh[i]=h->data[i];}
	for(i=h_len;i<plen;i++){	hh[i].re=0;	hh[i].im=0;}

	FFT(xx,fxx,plen);
	FFT(hh,fhh,plen);
	delete[] xx; xx=NULL;	delete[] hh; hh=NULL;

	for(i=0;i<plen;i++)
	{	
		fhh[i]=Conj(fhh[i]);
		fxx[i]=cMult(fxx[i],fhh[i]);
	}

	IFFT(fxx,fhh,plen);
	//选择fhh的最大值并返回对应下标
	float max_val=0;	int index=0;
	float tt=0;  Complex ctt;
	for(i=0;i<plen;i++)
	{	
		ctt=Conj(fhh[i]);
		y[i]=ctt;//输出的相关序列
		tt=Modulu(ctt);
		if(max_val<tt)	{	max_val=tt;		index=i;}
	}
	delete[] fxx; fxx=NULL;	delete[] fhh; fhh=NULL;

	float left,right;
	if(index==0) {	left=Modulu(y[plen-1]);		right=Modulu(y[index+1]);}
	else{	
		if(index==(plen-1)) {	left=Modulu(y[index-1]);	right=Modulu(y[0]);}
		else {	left=Modulu(y[index-1]);	right=Modulu(y[index+1]);}
	}
	float fa=(left+right-2*max_val)/2;
	float fb=(right-left)/2;
	double xstart=(-1)*fb/(2*fa);
	double shift=index+xstart;
	return shift;
}

///////////////
//矩阵基本运算
///////////////
//矩阵pC=pA*pB
//pC应该在主调函数中声明和分配空间
void Matrix_Mult(const Complex_Mat* pA,const Complex_Mat* pB,Complex_Mat* pC)
{
	assert((pA->mcol==pB->nrow) && (pA->nrow==pC->nrow) && (pB->mcol==pC->mcol));
	int i,j,k;
	Complex ss,tt;
	for(i=0;i<pC->nrow;i++)
	{
		for(j=0;j<pC->mcol;j++)
		{
			ss.re=0; ss.im=0;
			for(k=0;k<pA->mcol;k++)
			{
				tt=cMult(pA->data[i][k],pB->data[k][j]);
				ss=cAdd(ss,tt);
			}
			pC->data[i][j]=ss;
		}
	}
}

//pA矩阵的共轭转置(Conjugate Transpose).pB=pA'
//pB应该在主调函数中声明和分配空间
void Matrix_cTrans(const Complex_Mat* pA,Complex_Mat* pB)
{
	assert((pA->nrow==pB->mcol) && (pA->mcol==pB->nrow));
	for(int i=0;i<pB->nrow;i++)
	{
		for(int j=0;j<pB->mcol;j++)
		{	
			pB->data[i][j]=Conj(pA->data[j][i]);
		}
	}
}


//pA矩阵的转置(Transpose).pB=pA.'
//pB应该在主调函数中声明和分配空间
void Matrix_Trans(const Complex_Mat* pA,Complex_Mat* pB)
{
	assert((pA->nrow==pB->mcol) && (pA->mcol==pB->nrow));
	for(int i=0;i<pB->nrow;i++)
	{
		for(int j=0;j<pB->mcol;j++)
		{	
			pB->data[i][j]=pA->data[j][i];
		}
	}
}

//将矩阵pA以列优先转换为列矢量(即pB->mcol=1)
//pB应该在主调函数中声明和分配空间
void Mat2Vec(const Complex_Mat* pA,Complex_Mat* pB)
{
	int dim=(pA->nrow)*(pA->mcol);
	int st=0;
	assert((dim==pB->nrow) && (pB->mcol==1));
	for(int j=0;j<pA->mcol;j++)
	{
		st=j*(pA->nrow);
		for(int i=0;i<pA->nrow;i++)
		{	
			pB->data[st+i][0]=pA->data[i][j];
		}
	}
}

//对方阵pRx[nrow][nrow]求逆并存入inv_Rx
void Matrix_Inv(Complex** pRx,Complex** inv_Rx,int nrow)
{
	int mcol=nrow;//要求输入为方阵
	int n=nrow;//固定的迭代次数等于矩阵维数
	int psize=nrow*mcol;
	double* ar=new double[psize]; assert(ar!=NULL);//保存输入矩阵的实部(行优先).内部格式
	double* ai=new double[psize]; assert(ai!=NULL);//保存输入矩阵的虚部(行优先).内部格式

	int *is,*js,i,j,k,l,u,v,w;
    double p,q,s,t,d,b;
    is=new int[n]; assert(is!=NULL);//存放模值最大值元素的行标
    js=new int[n]; assert(js!=NULL);//存放模值最大值元素的列标
	
	for(i=0;i<n;i++) { is[i]=0;	js[i]=0;}

	for(i=0;i<nrow;i++)
	{
		for(j=0;j<mcol;j++)
		{
			u=i*nrow+j;//[i][j]
			ar[u]=pRx[i][j].re;
			ai[u]=pRx[i][j].im;
		}
	}//转换为内部格式
	//求逆过程
	for(k=0; k<n; k++)
	{
		d=0.0;
        for(i=k; i<n; i++)
		{
			for(j=k; j<n; j++)
			{ 
				u=i*n+j;//[i][j]
				p=ar[u]*ar[u]+ai[u]*ai[u];//模值
				if (p>d) { d=p; is[k]=i; js[k]=j;}
			}
		}
		if( (d<EPS) && (d>(-1)*EPS))
		{ 
			delete[] is; is=NULL;
			delete[] js; js=NULL;
			cout<<"err**not inv"<<endl;		exit(1);
		}
        if(is[k]!=k)
		{
          for(j=0; j<n; j++)
            { 
			  u=k*n+j; v=is[k]*n+j;
              t=ar[u]; ar[u]=ar[v]; ar[v]=t;
              t=ai[u]; ai[u]=ai[v]; ai[v]=t;
            }
		}
        if(js[k]!=k)
		{
          for(i=0; i<n; i++)
            { u=i*n+k; v=i*n+js[k];
              t=ar[u]; ar[u]=ar[v]; ar[v]=t;
              t=ai[u]; ai[u]=ai[v]; ai[v]=t;
            }
		}
        l=k*n+k;
        ar[l]=ar[l]/d; ai[l]=-ai[l]/d;
        for(j=0; j<=n-1; j++)
		{
          if(j!=k)
            { 
			  u=k*n+j;
              p=ar[u]*ar[l]; q=ai[u]*ai[l];
              s=(ar[u]+ai[u])*(ar[l]+ai[l]);
              ar[u]=p-q; ai[u]=s-p-q;
            }
		}
        for(i=0; i<n; i++)
		{
          if(i!=k)
		  {
			  v=i*n+k;
              for(j=0; j<n; j++)
			  {
                if(j!=k)
                  { 
					u=k*n+j;  w=i*n+j;
                    p=ar[u]*ar[v]; q=ai[u]*ai[v];
                    s=(ar[u]+ai[u])*(ar[v]+ai[v]);
                    t=p-q; b=s-p-q;
                    ar[w]=ar[w]-t;
                    ai[w]=ai[w]-b;
                  }
			  }
		  }
		}
        for(i=0; i<n; i++)
		{
          if(i!=k)
          { 
			  u=i*n+k;
              p=ar[u]*ar[l]; q=ai[u]*ai[l];
              s=(ar[u]+ai[u])*(ar[l]+ai[l]);
              ar[u]=q-p; ai[u]=p+q-s;
		  }
		}
	}

    for(k=n-1; k>=0; k--)
	{ 
		if(js[k]!=k)
		{
          for(j=0; j<=n-1; j++)
		  { 
			  u=k*n+j; v=js[k]*n+j;
              t=ar[u]; ar[u]=ar[v]; ar[v]=t;
              t=ai[u]; ai[u]=ai[v]; ai[v]=t;
		  }
		}
        if(is[k]!=k)
		{
          for(i=0; i<n; i++)
		  { 
			  u=i*n+k; v=i*n+is[k];
              t=ar[u]; ar[u]=ar[v]; ar[v]=t;
              t=ai[u]; ai[u]=ai[v]; ai[v]=t;
		  }
		}
	}//求逆完成

	//向输出inv_Rx赋值
	for(i=0;i<nrow;i++)
	{
		for(j=0;j<mcol;j++)
		{
			u=i*nrow+j;
			inv_Rx[i][j].re=(float) ar[u];
			inv_Rx[i][j].im=(float) ai[u];
		}
	}
	delete[] ar;	ar=NULL;	delete[] ai;	ai=NULL;
}

//重载方阵求逆Matrix_Inv
void Matrix_Inv(const Complex_Mat* pA,Complex_Mat* inv_A)
{
	assert( (pA->nrow==pA->mcol) && (pA->nrow==inv_A->nrow) && (pA->mcol==inv_A->mcol));
	Matrix_Inv(pA->data,inv_A->data,pA->nrow);
}
void Matrix_Inv(const dComplex_Mat* pA,dComplex_Mat* inv_A)
{
	assert( (pA->nrow==pA->mcol) && (pA->nrow==inv_A->nrow) && (pA->mcol==inv_A->mcol));
	Matrix_Inv(pA->data,inv_A->data,pA->nrow);
}
//求矩阵pRx[nrow][mcol]的伪逆并存入pinv_Rx
//如果nrow>mcol=>pinv_Rx=inv(pRx'*pRx)*pRx'==>pinv_Rx[mcol][nrow]==>pinv_Rx*pRx=I;
//如果nrow<mcol=>pinv_Rx=pRx'*inv(pRx*pRx')==>pinv_Rx[mcol][nrow]==>pRx*pinv_Rx=I;
//满足:dim(I)=min(nrow,mcol)
void Matrix_pInv(Complex** pRx,Complex** pinv_Rx,int nrow,int mcol)
{
	int dim=nrow<mcol?nrow:mcol;
	int psize=dim*dim;
	int i,j,k;
	Complex ss,tt;

	Complex** Rxx=new Complex*[psize]; assert(Rxx!=NULL);//pRx'*pRx or pRx*pRx'
	Init_pArray(Rxx,psize,psize);

	Complex** inv_Rxx=new Complex*[psize]; assert(inv_Rxx!=NULL);//
	Init_pArray(inv_Rxx,psize,psize);

	Complex** ct_Rx=new Complex*[mcol]; assert(ct_Rx!=NULL);//pRx的共轭转置距阵
	for(i=0;i<mcol;i++)
	{	
		ct_Rx[i]=new Complex[nrow]; assert(ct_Rx[i]!=NULL);
		for(j=0;j<nrow;j++){	ct_Rx[i][j]=Conj(pRx[j][i]);}//共轭转置
	}	

	if(nrow>mcol)
	{
		for(i=0;i<dim;i++)
		{
			for(j=0;j<dim;j++)
			{
				ss.re=0;	ss.im=0;
				for(k=0;k<nrow;k++)
				{
					tt=cMult(ct_Rx[i][k],pRx[k][j]);
					ss=cAdd(ss,tt);
				}//矩阵乘法
				Rxx[i][j]=ss;
			}
		}
		Matrix_Inv(Rxx,inv_Rxx,dim);	
		for(i=0;i<dim;i++)
		{
			for(j=0;j<nrow;j++)
			{
				ss.re=0;	ss.im=0;
				for(k=0;k<dim;k++)
				{
					tt=cMult(inv_Rxx[i][k],ct_Rx[k][j]);
					ss=cAdd(ss,tt);
				}//矩阵乘法
				pinv_Rx[i][j]=ss;
			}
		}		
		Discard_pArray(Rxx,dim);	delete[] Rxx;	Rxx=NULL;
		Discard_pArray(inv_Rxx,dim);	delete[] inv_Rxx; inv_Rxx=NULL;
		Discard_pArray(ct_Rx,mcol);	delete[] ct_Rx; ct_Rx=NULL;
	}
	else
	{
		if(nrow<mcol)
		{
			for(i=0;i<dim;i++)
			{
				for(j=0;j<dim;j++)
				{
					ss.re=0;	ss.im=0;
					for(k=0;k<mcol;k++)
					{
						tt=cMult(pRx[i][k],ct_Rx[k][j]);
						ss=cAdd(ss,tt);
					}//矩阵乘法
					Rxx[i][j]=ss;
				}
			}
			Matrix_Inv(Rxx,inv_Rxx,dim);
			//inv_Rxx*pRx
			for(i=0;i<mcol;i++)
			{
				for(j=0;j<dim;j++)
				{
					ss.re=0;	ss.im=0;
					for(k=0;k<dim;k++)
					{
						tt=cMult(ct_Rx[i][k],inv_Rxx[k][j]);
						ss=cAdd(ss,tt);
					}//矩阵乘法
					pinv_Rx[i][j]=ss;
				}
			}
			Discard_pArray(Rxx,dim);	delete[] Rxx;	Rxx=NULL;
			Discard_pArray(inv_Rxx,dim);	delete[] inv_Rxx; inv_Rxx=NULL;
			Discard_pArray(ct_Rx,mcol);	delete[] ct_Rx; ct_Rx=NULL;
		}
		else
		{  
			Matrix_Inv(pRx,pinv_Rx,dim);
			Discard_pArray(Rxx,dim);	delete[] Rxx;	Rxx=NULL;
			Discard_pArray(inv_Rxx,dim);	delete[] inv_Rxx; inv_Rxx=NULL;
			Discard_pArray(ct_Rx,mcol);	delete[] ct_Rx; ct_Rx=NULL;
		}
	}
}
//重载矩阵伪逆Matrix_pInv
void Matrix_pInv(const Complex_Mat* pA,Complex_Mat* pinv_A)
{
	assert((pinv_A->mcol==pA->nrow) && (pinv_A->nrow==pA->mcol));
	Matrix_pInv(pA->data,pinv_A->data,pA->nrow,pA->mcol);
}

//////////////////////
//基本统计量估计函数
//////////////////////
//返回数据的均值
double mean(const sqList* sq)
{
	double mm=0;
	int num=0;
	sqList_node* q=sq->head->next;
	while(q!=NULL)
	{	
		mm+=q->val;
		q=q->next;
		num++;
	}
	mm=mm/num;
	return mm;
}
//返回n点双精度数组的均值
double mean(const double* pd,int n)
{
	double mm=0;
	for(int i=0;i<n;i++) {	mm+=pd[i];}
	mm=mm/n;
	return mm;
}
//返回n点双精度数组的均值
float mean(const float* pd,int n)
{
	float mm=0;
	for(int i=0;i<n;i++) {	mm+=pd[i];}
	mm=mm/n;
	return mm;
}
//返回n点复数数组的均值
Complex mean(const Complex* pd,int n)
{
	Complex mm;
	mm.re=0;	mm.im=0;
	for(int i=0;i<n;i++) {	mm.re+=pd[i].re;	mm.im+=pd[i].im;}
	mm.re=mm.re/n;	mm.im=mm.im/n;
	return mm;
}
	
//返回标准差
double std(const sqList* sq)
{
	double mm=mean(sq);
	double sum=0;
	double tt=0;
	int num=0;
	sqList_node* q=sq->head->next;
	while(q!=NULL)
	{	
		tt=q->val-mm;
		sum+=(tt*tt);
		q=q->next;
		num++;
	}
	sum=sqrt(sum/(num-1));
	return sum;
}

//返回中值
double median(const sqList* sq)
{
	double* dp=new double[sq->len];//值数组
	int* index=new int[sq->len];//下标数组
	sqList_node* q=sq->head->next;
	int i=0;
	while(q!=NULL)
	{
		dp[i]=q->val;	index[i]=q->index;
		q=q->next;
		i++;
	}
	//排序
	double td;
	int tindex,j,low,high,mid;
	for(i=1;i<sq->len;i++)
	{
		td=dp[i];	tindex=index[i];
		low=0; high=i-1; mid=(low+high)/2;
		while(low<=high)
		{	
			if(dp[mid]>td)	{	high=mid-1;}
			else	{	low=mid+1;}
			mid=(low+high)/2;
		}//查找插入位置
		for(j=i-1;j>=low;j--)
		{
			dp[j+1]=dp[j];
			index[j+1]=index[j];
		}
		dp[low]=td;
		index[low]=tindex;		
	}
	tindex=ceil(sq->len,2)-1;
	td=dp[tindex];

	delete[] dp; dp=NULL;	delete[] index; index=NULL;
	return td;
}
//返回n点数组的中值
double median(const double* pd,int n)
{
	double* dp=new double[n];//值数组
	int* index=new int[n];//下标数组
	int i=0;
	for(i=0;i<n;i++) {	dp[i]=pd[i];	index[i]=i;}
	//排序
	double td;
	int tindex,j,low,high,mid;
	for(i=1;i<n;i++)
	{
		td=dp[i];	tindex=index[i];
		low=0; high=i-1; mid=(low+high)/2;
		while(low<=high)
		{	
			if(dp[mid]>td)	{	high=mid-1;}
			else	{	low=mid+1;}
			mid=(low+high)/2;
		}//查找插入位置
		for(j=i-1;j>=low;j--)
		{
			dp[j+1]=dp[j];
			index[j+1]=index[j];
		}
		dp[low]=td;
		index[low]=tindex;		
	}
	tindex=ceil(n,2)-1;
	td=dp[tindex];
	delete[] dp; dp=NULL;
	delete[] index; index=NULL;
	return td;
}

//对n点数组pd的元素从小到大升序排序
//返回排序后的数组ps和对应元素在pd中的下标
void sort(const double* pd,double* ps,int* index,int n)
{
	int i,j,low,high,mid,tindex;
	for(i=0;i<n;i++) {	ps[i]=pd[i];	index[i]=i;}
	//排序
	double td;
	for(i=1;i<n;i++)
	{
		td=ps[i];	tindex=index[i];
		low=0; high=i-1; mid=(low+high)/2;
		while(low<=high)
		{	
			if(ps[mid]>td)	{	high=mid-1;}
			else	{	low=mid+1;}
			mid=(low+high)/2;
		}//查找插入位置
		for(j=i-1;j>=low;j--)
		{
			ps[j+1]=ps[j];
			index[j+1]=index[j];
		}
		ps[low]=td;
		index[low]=tindex;		
	}
}


//返回矩阵元素模值的最大值和对应下标
//index[0]存放行下标,index[1]存放列下标
float Find_Max_Value(const Complex_Mat* in_ptr,int* index)
{
	float  max_val,temp;
	max_val=Modulu(in_ptr->data[0][0]); index[0]=0; index[1]=0;
	for(int i=0;i<in_ptr->nrow;i++)
	{
		for(int j=0;j<in_ptr->mcol;j++)
		{
			temp=Modulu(in_ptr->data[i][j]);
			if(max_val < temp)
			{	index[0]=i;	index[1]=j;	max_val=temp;}
		}
	}
	return max_val;
}
//返回矩阵元素模值的最大值和对应下标(存放在index中)
float Find_Max_Value(const Complex_Vec* in_ptr,int& index)
{
	float  max_val,temp;
	max_val=Modulu(in_ptr->data[0]); index=0;
	for(int i=0;i<in_ptr->len;i++)
	{
		temp=Modulu(in_ptr->data[i]);
		if(max_val < temp)
		{	index=i;	max_val=temp;}
	}
	return max_val;
}