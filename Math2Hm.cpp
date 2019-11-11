#include "CFAR_SubFunction.h"


/////////////////////////////////////////////////
//ʵ��GMTI���ݷ�����ڲ��������ݸ�ʽ
/////////////////////////////////////////////////
//   $zwYang, The RSP.
//   $Revision: 1.0 $  $Date: 2009/04/14 18:20 $
/////////////////////////////////////////////////
#define eb_eps 1e-5 //���߹��Ƶ��ڲ�����
#define TRUE 1
#define FALSE 0

////////////////////////////
//GMTI���ݴ����ڲ���������
///////////////////////////
//����
Complex::Complex() { re=0; im=0;}
Complex::Complex(float re,float im)
{ this->re=re; this->im=im;}
void Complex::Show()
{ cout<<re<<"+i"<<im<<endl;}


//��ʸ��
Complex_Vec::Complex_Vec(int max_size)
{
	if(max_size > 0)
	{
		data=new Complex[max_size];	assert(data!=NULL);//��������������ֹ
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

//������
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
}//����һ�����.data[nrow][mcol]
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
}//���췽��

Complex_Mat::~Complex_Mat()
{
	for(int n=0;n<nrow;n++)
	{
		delete[] data[n];
		data[n]=NULL;
	}
	delete[] data; data=NULL;
}//��������

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
//ʵ�ֹ����������źŴ����Ӻ���
/////////////////////////////////////////////////
//   $zwYang, The RSP.
//   $Revision: 1.0 $  $Date: 2009/04/16 17:20 $
/////////////////////////////////////////////////

////////////////
//ʵ����������
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

//��ָ��ת�ɸ���
extern inline Complex CompExp2Comp(Complex_Exp x)
{
	Complex y;
	y.re=(float) (x.rho*cos(x.theta));
	y.im=(float) (x.rho*sin(x.theta));
	return y;
}

//���ظ�����ģֵ
extern inline float Modulu(Complex x)
{
	float r=(float) (sqrt(x.re*x.re+x.im*x.im));
	return r;
}

//���ظ����Ĺ���ֵ
extern inline Complex Conj(Complex x)
{
	Complex y;
	y.re=x.re; y.im=(-1)*x.im;
	return y;
}

//��x/y�Ľ������ȡ��.ceil(5,4)=2
extern inline int ceil(int x,int y)
{
	int s=x/y;
	if((s*y) < x){	s=s+1;}
	return s;
}

//��x/y�Ľ������ȡ��.fix(5,2)=2
extern inline int fix(int x,int y)
{
	int s=x/y;
	if((s*y) > x){	s=s-1;}
	return s;
}

//////////////////////
//ָ������������ͷ� 
//////////////////////
//�������n��Ԫ�ص�ָ��������ָ�������ÿ��Ԫ��ָ�����m��Ԫ�ص�һά����
//pArray=new Complex* [n]����Ӧ������������������
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

//�ͷ�n��Ԫ�ص�ָ������pArray��ÿһ��Ԫ��(��ӦΪÿ��һά����)
//pArray������ͷ�Ӧ�����������������
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
//����ת��Ϊ��ָ����ȡ��ǲ���(-PI-PI)
////////////////////////////////////////
//����ת�ɸ�ָ��
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

//���ظ��������(rad.[-pi:pi))
float Angle(Complex x)
{
	
	double theta; 
	int sign=1;//x.im�ķ���.1:��ʾ��;-1��ʾ��
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
//2�Ķ�������ȡ��
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

//���Ϊ����Ķ�����������ֵ����bitsλ����
//���磺Bit_Inv(6,5):6->00110->01100=12
int Bit_Inv(int x,int bits)
{
	int s=0; //sΪҪ���ص��� 
	int index=bits-1; 
	while(index>=0) 
	{   
		s+=(x%2)*(int)pow(2,index); 
		x=x/2; 
		index--;
	} 
	return s; 
}//

//����������Ԫ�ذ��±���ж������������ŵõ��µ����� 
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


//���ٸ���Ҷ�任(FFT)
//��len����������н���plen��FFT����,�����out_ptrָ���ҳ���Ϊplen
//���plen<len����������н��н�ȡ�����plen>len����������������ӳ�
//plen Ϊ2����������(plen=2^n)����ÿ����㷨��������DFT����
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
	}//���ݱ任����plen��ԭ���н��н�ȡ����

	FFT(inn_iptr,out_ptr,plen);//����plen��FFT
    delete[] inn_iptr; inn_iptr=NULL;
}

//���ٸ���Ҷ�任
//��len����������н���len��FFT����,�����out_ptrָ����������ȳ�
//lenΪ2����������(len=2^n)�����ʱ���2��ȡ�����㷨��������DFT����
void FFT(const Complex* in_ptr, Complex* out_ptr, int len)
{
	int stage=log2(len);
	int flag=(int)pow(2,stage);
	if( flag==len){	flag=1;}//len=2^n
	else {flag=0;}

	double pha=((-2)*PI/len);//DFT��Ƶ�ʼ��(rad)
	Complex ss,dd,tt;//

	if(flag==1)
	{
		int pn=len/2;//ÿһ���ĵ����ܵ�Ԫ��
		Complex* W_f=new Complex[pn]; assert(W_f!=NULL);//��ת����
		for(int i=0;i<pn;i++)
		{
			W_f[i].re=(float)cos(pha*i);
			W_f[i].im=(float)sin(pha*i);
		}//������ת����

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
	}//len=2^n���ÿ����㷨
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
	}//DFT���㷽ʽ
}//reload FFT function


//��in_ptr����άFFT���ѽ�����浽out_ptr
void FFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr)
{
	assert((in_ptr->nrow==out_ptr->nrow)&&(in_ptr->mcol==out_ptr->mcol));	
	Complex* xx=new Complex[in_ptr->nrow]; assert(xx!=NULL);//������
	Complex* fxx=new Complex[in_ptr->nrow]; assert(fxx!=NULL);//�任����

	for (int i=0;i<in_ptr->nrow;i++)
	{	FFT(in_ptr->data[i],out_ptr->data[i],in_ptr->mcol);}//���б任

	for (i=0;i<out_ptr->mcol;i++)
	{
		for (int j=0;j<out_ptr->nrow;j++){	xx[j]=out_ptr->data[j][i];}
		FFT(xx,fxx,in_ptr->nrow);
		for (j=0;j<out_ptr->nrow;j++){	out_ptr->data[j][i]=fxx[j];}
	}//���б任
	delete[] xx; xx=NULL;		delete[] fxx; fxx=NULL;
}

//�����渵��Ҷ�任(IFFT)
//��len����������н���len��IFFT����,�����out_ptrָ����������ȳ�
//lenΪ2����������(len=2^n)�����ʱ���2��ȡ�����㷨��������IDFT����
void IFFT(const Complex* in_ptr, Complex* out_ptr,int len)
{
	int stage=log2(len);
	int flag=(int)pow(2,stage);
	if( flag==len){	flag=1;}//len=2^n
	else {flag=0;}

	double pha=((2)*PI/len);//IDFT��Ƶ�ʼ��(rad)
	Complex ss,dd,tt;//

	if(flag==1)
	{
		int pn=len/2;//ÿһ���ĵ����ܵ�Ԫ��
		Complex* W_f=new Complex[pn]; assert(W_f!=NULL);//��ת����
		for(int i=0;i<pn;i++)
		{
			W_f[i].re=(float)cos(pha*i);
			W_f[i].im=(float)sin(pha*i);
		}//������ת����

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
	}//len=2^n���ÿ����㷨
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
	}//IDFT���㷽ʽ
}


//��in_ptr����άIFFT���ѽ�����浽out_ptr
void IFFT2(const Complex_Mat* in_ptr, Complex_Mat* out_ptr)
{
	assert((in_ptr->nrow==out_ptr->nrow)&&(in_ptr->mcol==out_ptr->mcol));
	Complex* xx=new Complex[in_ptr->nrow]; assert(xx!=NULL);//������
	Complex* fxx=new Complex[in_ptr->nrow]; assert(fxx!=NULL);//�任����

	for (int i=0;i<in_ptr->nrow;i++)
	{ IFFT(in_ptr->data[i],out_ptr->data[i],in_ptr->mcol);}//���б任

	for (i=0;i<out_ptr->mcol;i++)
	{
		for (int j=0;j<out_ptr->nrow;j++){ xx[j]=out_ptr->data[j][i];}	
		IFFT(xx,fxx,in_ptr->nrow);
		for (j=0;j<out_ptr->nrow;j++){ out_ptr->data[j][i]=fxx[j];}
	}//���б任
	delete[] xx; xx=NULL;		delete[] fxx; fxx=NULL;
}

//��len����������in_ptr����ǰ�������ֺ����out_ptr
void FFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len)
{
	int pos=ceil(len,2);//ǰ�������ֵķֽ���±�
	int i,j=0;
	for(i=pos;i<len;i++){ out_ptr[j]=in_ptr[i];	j++;}//��벿��
	for(i=0;i<pos;i++){ out_ptr[j]=in_ptr[i]; j++;}//ǰ�벿��
}

//������in_ptr��ǰ��(������)�����ֽ�����������out_ptr��
void FFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr)
{
	assert(in_ptr->len==out_ptr->len);
	int pos=ceil(in_ptr->len,2);
	int i,j=0;
	for(i=pos;i<in_ptr->len;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
	for(i=0;i<pos;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
}

//�Ծ���in_ptr��ǰ�����������ֽ����������out_ptr��
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
	}//������(������)
	else if (in_ptr->mcol==1)
	{
		pos=ceil(in_ptr->nrow,2); j=0;
		for(i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
		for(i=0;i<pos;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
	}//������(������)
	else
	{
		pos=ceil(in_ptr->nrow,2);
		for(int p=0;p<in_ptr->mcol;p++)
		{
			j=0;
			for (i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][p]=in_ptr->data[i][p]; j++;}
			for (i=0;i<pos;i++){ out_ptr->data[j][p]=in_ptr->data[i][p];j++;}	
		}//ÿһ��ǰ�������ֽ�����Ľ������out_ptr��Ӧ��
		pos=ceil(in_ptr->mcol,2);
		for(p=0;p<in_ptr->nrow;p++)
		{
			j=0;
			for(i=pos;i<in_ptr->mcol;i++){ temp[j]=out_ptr->data[p][i]; j++;}
			for(i=0;i<pos;i++){ temp[j]=out_ptr->data[p][i]; j++;}//�ݴ���temp
			for(i=0;i<in_ptr->mcol;i++){ out_ptr->data[p][i]=temp[i];}//����out_ptr
		}//ÿһ��ǰ�������ֽ�����Ľ������out_ptr��Ӧ��
	}//
	delete[] temp; temp=NULL;
}

//��io_ptr��ָ�����������ǰ�����������ֽ����������������
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
		for(i=0;i<io_ptr->nrow;i++){ io_ptr->data[i][p]=temp[i];}//����io_ptr
	}//ÿһ��ǰ�������ֽ�����Ľ������io_ptr��Ӧ��
	pos=ceil(io_ptr->mcol,2);
	for(p=0;p<io_ptr->nrow;p++)
	{
		j=0;
		for(i=pos;i<io_ptr->mcol;i++){ temp[j]=io_ptr->data[p][i]; j++;}
		for(i=0;i<pos;i++){ temp[j]=io_ptr->data[p][i]; j++;}//�ݴ���temp
		for(i=0;i<io_ptr->mcol;i++){ io_ptr->data[p][i]=temp[i];}//����out_ptr
	}//ÿһ��ǰ�������ֽ�����Ľ������io_ptr��Ӧ��	
	delete[] temp; temp=NULL;
}

//��len����������in_ptr����ǰ�������ֺ����out_ptr
//IFFTSHIFT(FFTSHIFT(X))=X;
void IFFTSHIFT(const Complex* in_ptr, Complex* out_ptr, int len)
{
	int pos=fix(len,2);
	int i,j=0;
	for(i=pos;i<len;i++){ out_ptr[j]=in_ptr[i]; j++;}//��벿��
	for(i=0;i<pos;i++){ out_ptr[j]=in_ptr[i]; j++;}//ǰ�벿��
}

//����IFFTSHIFT
void IFFTSHIFT(const Complex_Vec* in_ptr, Complex_Vec* out_ptr)
{
	assert(in_ptr->len==out_ptr->len);
	int pos=fix(in_ptr->len,2);
	int i,j=0;
	for(i=pos;i<in_ptr->len;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
	for(i=0;i<pos;i++){ out_ptr->data[j]=in_ptr->data[i]; j++;}
}

//����IFFTSHIFT
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
	}//������(������)
	else if (in_ptr->mcol==1)
	{
		pos=fix(in_ptr->nrow,2); j=0;
		for(i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
		for(i=0;i<pos;i++){ out_ptr->data[j][0]=in_ptr->data[i][0]; j++;}
	}//������(������)
	else
	{
		pos=fix(in_ptr->nrow,2);
		for(int p=0;p<in_ptr->mcol;p++)
		{
			j=0;
			for (i=pos;i<in_ptr->nrow;i++){ out_ptr->data[j][p]=in_ptr->data[i][p]; j++;}
			for (i=0;i<pos;i++){ out_ptr->data[j][p]=in_ptr->data[i][p];j++;}	
		}//ÿһ��ǰ�������ֽ�����Ľ������out_ptr��Ӧ��
		pos=fix(in_ptr->mcol,2);
		for(p=0;p<in_ptr->nrow;p++)
		{
			j=0;
			for(i=pos;i<in_ptr->mcol;i++){ temp[j]=out_ptr->data[p][i]; j++;}
			for(i=0;i<pos;i++){ temp[j]=out_ptr->data[p][i]; j++;}//�ݴ���temp
			for(i=0;i<in_ptr->mcol;i++){ out_ptr->data[p][i]=temp[i];}//����out_ptr
		}//ÿһ��ǰ�������ֽ�����Ľ������out_ptr��Ӧ��
	}//
	delete[] temp; temp=NULL;
}

//��io_ptr��ָ�����������ǰ�����������ֽ����������������
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
		for(i=0;i<io_ptr->nrow;i++){ io_ptr->data[i][p]=temp[i];}//����io_ptr
	}//ÿһ��ǰ�������ֽ�����Ľ������io_ptr��Ӧ��
	pos=fix(io_ptr->mcol,2);
	for(p=0;p<io_ptr->nrow;p++)
	{
		j=0;
		for(i=pos;i<io_ptr->mcol;i++){ temp[j]=io_ptr->data[p][i]; j++;}
		for(i=0;i<pos;i++){ temp[j]=io_ptr->data[p][i]; j++;}//�ݴ���temp
		for(i=0;i<io_ptr->mcol;i++){ io_ptr->data[p][i]=temp[i];}//����out_ptr
	}//ÿһ��ǰ�������ֽ�����Ľ������io_ptr��Ӧ��	
	delete[] temp; temp=NULL;
}

//�������ξ��(����FFT��IFFT)
//��������x��h
//�������y=x*h(����Ϊlength(x)+length(h)-1)
//���������FFT��IFFT�������ξ��
void Conv(const Complex_Vec* x,const Complex_Vec*h,Complex* y)
{
	int x_len,h_len,y_len,plen;
	x_len=x->len;	h_len=h->len;	y_len=x_len+h_len-1;
	plen=log2(y_len);	plen=(int) pow(2,plen);//FFT�ı任����

	Complex* xx=new Complex[plen];  assert(xx!=NULL);//��������������
	Complex* hh=new Complex[plen];  assert(hh!=NULL);//��������������
	Complex* fxx=new Complex[plen];	assert(fxx!=NULL);//FFT�ı任���
	Complex* fhh=new Complex[plen];	assert(fhh!=NULL);//FFT�ı任���
	
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

	for(i=0;i<y_len;i++){	y[i]=fhh[i];}//��ȡ�������Чֵ

	delete[] fxx; fxx=NULL;	 delete[] fhh; fhh=NULL;
}


//������������x��h��ѭ�����(����FFT��IFFT)
//����h�����x��ƫ����(���ֵ��Ӧ��λ��)��ѭ���������y
//�������y(m)=x(n).conj(h(n-m))--h(n)����ѭ����λ
//      ��y(m)=conj(IFFT(FFT(x,len).*conj(FFT(h,len))))
//������г���lenΪx��h����󳤶�ֵ.m=[0:len-1]
//��������жԽ϶����в������ϳ����еĳ���
double Corr(const Complex_Vec* x,const Complex_Vec* h,Complex* y)
{
	int x_len=x->len;
	int h_len=h->len;
	int plen=(x_len>h_len?x_len:h_len);

	Complex* xx=new Complex[plen];  assert(xx!=NULL);//��������������
	Complex* hh=new Complex[plen];  assert(hh!=NULL);//��������������
	Complex* fxx=new Complex[plen];	assert(fxx!=NULL);//FFT�ı任���
	Complex* fhh=new Complex[plen];	assert(fhh!=NULL);//FFT�ı任���
	
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
	//ѡ��fhh�����ֵ�����ض�Ӧ�±�
	float max_val=0;	int index=0;
	float tt=0;  Complex ctt;
	for(i=0;i<plen;i++)
	{	
		ctt=Conj(fhh[i]);
		y[i]=ctt;//������������
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
//�����������
///////////////
//����pC=pA*pB
//pCӦ�������������������ͷ���ռ�
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

//pA����Ĺ���ת��(Conjugate Transpose).pB=pA'
//pBӦ�������������������ͷ���ռ�
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


//pA�����ת��(Transpose).pB=pA.'
//pBӦ�������������������ͷ���ռ�
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

//������pA��������ת��Ϊ��ʸ��(��pB->mcol=1)
//pBӦ�������������������ͷ���ռ�
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

//�Է���pRx[nrow][nrow]���沢����inv_Rx
void Matrix_Inv(Complex** pRx,Complex** inv_Rx,int nrow)
{
	int mcol=nrow;//Ҫ������Ϊ����
	int n=nrow;//�̶��ĵ����������ھ���ά��
	int psize=nrow*mcol;
	double* ar=new double[psize]; assert(ar!=NULL);//������������ʵ��(������).�ڲ���ʽ
	double* ai=new double[psize]; assert(ai!=NULL);//�������������鲿(������).�ڲ���ʽ

	int *is,*js,i,j,k,l,u,v,w;
    double p,q,s,t,d,b;
    is=new int[n]; assert(is!=NULL);//���ģֵ���ֵԪ�ص��б�
    js=new int[n]; assert(js!=NULL);//���ģֵ���ֵԪ�ص��б�
	
	for(i=0;i<n;i++) { is[i]=0;	js[i]=0;}

	for(i=0;i<nrow;i++)
	{
		for(j=0;j<mcol;j++)
		{
			u=i*nrow+j;//[i][j]
			ar[u]=pRx[i][j].re;
			ai[u]=pRx[i][j].im;
		}
	}//ת��Ϊ�ڲ���ʽ
	//�������
	for(k=0; k<n; k++)
	{
		d=0.0;
        for(i=k; i<n; i++)
		{
			for(j=k; j<n; j++)
			{ 
				u=i*n+j;//[i][j]
				p=ar[u]*ar[u]+ai[u]*ai[u];//ģֵ
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
	}//�������

	//�����inv_Rx��ֵ
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

//���ط�������Matrix_Inv
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
//�����pRx[nrow][mcol]��α�沢����pinv_Rx
//���nrow>mcol=>pinv_Rx=inv(pRx'*pRx)*pRx'==>pinv_Rx[mcol][nrow]==>pinv_Rx*pRx=I;
//���nrow<mcol=>pinv_Rx=pRx'*inv(pRx*pRx')==>pinv_Rx[mcol][nrow]==>pRx*pinv_Rx=I;
//����:dim(I)=min(nrow,mcol)
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

	Complex** ct_Rx=new Complex*[mcol]; assert(ct_Rx!=NULL);//pRx�Ĺ���ת�þ���
	for(i=0;i<mcol;i++)
	{	
		ct_Rx[i]=new Complex[nrow]; assert(ct_Rx[i]!=NULL);
		for(j=0;j<nrow;j++){	ct_Rx[i][j]=Conj(pRx[j][i]);}//����ת��
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
				}//����˷�
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
				}//����˷�
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
					}//����˷�
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
					}//����˷�
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
//���ؾ���α��Matrix_pInv
void Matrix_pInv(const Complex_Mat* pA,Complex_Mat* pinv_A)
{
	assert((pinv_A->mcol==pA->nrow) && (pinv_A->nrow==pA->mcol));
	Matrix_pInv(pA->data,pinv_A->data,pA->nrow,pA->mcol);
}

//////////////////////
//����ͳ�������ƺ���
//////////////////////
//�������ݵľ�ֵ
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
//����n��˫��������ľ�ֵ
double mean(const double* pd,int n)
{
	double mm=0;
	for(int i=0;i<n;i++) {	mm+=pd[i];}
	mm=mm/n;
	return mm;
}
//����n��˫��������ľ�ֵ
float mean(const float* pd,int n)
{
	float mm=0;
	for(int i=0;i<n;i++) {	mm+=pd[i];}
	mm=mm/n;
	return mm;
}
//����n�㸴������ľ�ֵ
Complex mean(const Complex* pd,int n)
{
	Complex mm;
	mm.re=0;	mm.im=0;
	for(int i=0;i<n;i++) {	mm.re+=pd[i].re;	mm.im+=pd[i].im;}
	mm.re=mm.re/n;	mm.im=mm.im/n;
	return mm;
}
	
//���ر�׼��
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

//������ֵ
double median(const sqList* sq)
{
	double* dp=new double[sq->len];//ֵ����
	int* index=new int[sq->len];//�±�����
	sqList_node* q=sq->head->next;
	int i=0;
	while(q!=NULL)
	{
		dp[i]=q->val;	index[i]=q->index;
		q=q->next;
		i++;
	}
	//����
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
		}//���Ҳ���λ��
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
//����n���������ֵ
double median(const double* pd,int n)
{
	double* dp=new double[n];//ֵ����
	int* index=new int[n];//�±�����
	int i=0;
	for(i=0;i<n;i++) {	dp[i]=pd[i];	index[i]=i;}
	//����
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
		}//���Ҳ���λ��
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

//��n������pd��Ԫ�ش�С������������
//��������������ps�Ͷ�ӦԪ����pd�е��±�
void sort(const double* pd,double* ps,int* index,int n)
{
	int i,j,low,high,mid,tindex;
	for(i=0;i<n;i++) {	ps[i]=pd[i];	index[i]=i;}
	//����
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
		}//���Ҳ���λ��
		for(j=i-1;j>=low;j--)
		{
			ps[j+1]=ps[j];
			index[j+1]=index[j];
		}
		ps[low]=td;
		index[low]=tindex;		
	}
}


//���ؾ���Ԫ��ģֵ�����ֵ�Ͷ�Ӧ�±�
//index[0]������±�,index[1]������±�
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
//���ؾ���Ԫ��ģֵ�����ֵ�Ͷ�Ӧ�±�(�����index��)
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