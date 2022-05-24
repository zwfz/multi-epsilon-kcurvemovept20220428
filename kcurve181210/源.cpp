#include<gl/glut.h>
#include<vector>
#include <Eigen/LU>
#include<Eigen/Dense>
#include <Eigen/StdVector>
#include<iostream>
#include<math.h>
#include<fstream>
#include <chrono>		//时间库

#include "tinyfiledialogs.h"
#include "utils1.h"

//#include <glew.h>
//#include <freeglut.h>
//#include <freeglut_ext.h>
/*使用二分法求ti 20220513*/
/*实现画曲率直方图 20220514*/
/*改用牛顿法求ti 20220523*/


using namespace std;
using namespace Eigen;

typedef Eigen::Vector2d EVec2d;
typedef vector<EVec2d, Eigen::aligned_allocator<EVec2d>> vecEg2dd;

typedef vector<vecEg2dd, Eigen::aligned_allocator<vecEg2dd>> VvecEg2dd;			//二维点列的矩阵容器
typedef vector<vector<double>> Vvector;			//double的矩阵容器

#define PI 3.14159265;

//vector<EVec2d> cps;
//vector<EVec2d> ci0, ci1, ci2;
vecEg2dd cps,ci0,ci1,ci2,cti,dcti;

vecEg2dd dcps;

VvecEg2dd MtrxCps_o,MtrxCps_c;				//控制点矩阵

vector<double> ti;
vector<double> lambdai;
vector<double> theta;
vector<double> ai;

Vvector Mtrxai_o,Mtrxai_c;				//ai的矩阵

double defaultai = 2. / 3.;

GLsizei winwidth = 800;
GLsizei winHeight = 800;
bool mouseLeftDown;
bool mouseRightDown;
int Num,i0;
bool flag0 = false;
bool closedflag = false;

bool aiflag = false;

void closed_ekcurve(vecEg2dd& cpsi, vector<double>& ai);
void opened_ekcurve(vecEg2dd& cpsi, vector<double>& ai);

void calculate_ci02(int n);
void opencalculate_ci02(int n);
void calculate_ci1(int n, vecEg2dd& cpsi, vector<double>& ai);
void opencalculate_ci1(int n, vecEg2dd& cpsi, vector<double>& ai);
void calculate_ti(int n, vecEg2dd& cpsi, vector<double>& ai);
void opencalculate_ti(int n, vecEg2dd& cpsi, vector<double>& ai);
void calculate_lambdai(int n, vector<double>& ai);
void opencalculate_lambdai(int n, vector<double>& ai);
double TrgArea(const EVec2d&p0, const EVec2d&p1, const EVec2d&p2);

bool closedCheckCvrg(int n, vecEg2dd& cpsi, vector<double>& ai);
bool openCheckCvrg(int n, vecEg2dd& cpsi, vector<double>& ai);

bool CheckItr(const vecEg2dd &prevci1,int n);
void drawQuBzr(EVec2d &ci0, EVec2d &ci1, EVec2d &ci2, double ai);
void drawCrvtrVctrClmn(EVec2d& ci0, EVec2d& ci1, EVec2d& ci2, double ai);

//double bisection(int i, EVec2d& cpsii, double aii);
double dCrvtr(int i, double ti, EVec2d& cpsii, double aii);
double ddCrvtr(int i, double ti, EVec2d& cpsii, double aii);
double Newton(int i, EVec2d& cpsii, double aii);

bool checkpt(int x, int y);
//void mouseMotionPT(int x, int y);
void OnMouse(int button, int state, int x, int y);
void myKeyboard(unsigned char key, int x, int y);
void MouseMove(int x, int y);

void menuFunc(int value);
void saveandload(int value);

void InpPtsFl();
void SvPtsFl();

void changeai(int i);
int getaiIndex(int x, int y);

void endcurve();

void myInit()
{
	
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.0, 0.0, 0.0, 0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0,800,0, 800);
}

void myReshape(int w,int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);
}

void endcurve()
{
	if (closedflag)
	{
		MtrxCps_c.push_back({ cps });			//闭式eK-curve点矩阵
		Mtrxai_c.push_back({ ai });				//ai值矩阵
	}
	if (!closedflag)
	{
		MtrxCps_o.push_back({ cps });			//开式eK-curve点矩阵
		Mtrxai_o.push_back({ ai });				//ai值矩阵
	}

	flag0 = false;
	cps.clear();
	ai.clear();
}

void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	Num = (int)cps.size();

	/*for (EVec2d& point : cps)
	{
		cout << "cps:" << point.transpose() << endl;
	}*/

	glColor3f(1, 1, 0);		//黄色
	glPointSize(7);
	glBegin(GL_POINTS);

	for (int i = 0; i < Num; ++i)
	{
		glVertex2d(cps[i][0], cps[i][1]);
	}

	glEnd();

	/*if (Num > 0)
	{
		Inputai();
	}*/
	ai.resize(Num, defaultai);


	glColor3f(0, 1, 1);		//蓝色
	glLineWidth(3);
	glBegin(GL_LINE_STRIP);   //画线
	for (int i = 0; i < Num; ++i)
	{
		glVertex2d(cps[i][0], cps[i][1]);
	}
	glEnd();

	
	if (Num > 2)
	{
		if ((Num>3) && closedflag)
		{
			closed_ekcurve(cps,ai);
			for (int j = 0; j < Num; ++j)
			{
				drawQuBzr(ci0[j], ci1[j], ci2[j],ai[j]);
				drawCrvtrVctrClmn(ci0[j], ci1[j], ci2[j], ai[j]);
			}
		}
		if (!closedflag)
		{
			opened_ekcurve(cps,ai);
			for (int j = 0; j < Num-2; ++j)
			{
				drawQuBzr(ci0[j], ci1[j], ci2[j],ai[j]);
				drawCrvtrVctrClmn(ci0[j], ci1[j], ci2[j], ai[j]);
			}
		}
		

		
	}

	//屏幕显示ai值
	glColor3d(1., 0.502, 0);		//橘黄
	if ((Num > 3) && closedflag)				//闭式ek-curve
	{

		for (int i = 0; i < Num; ++i)
		{
			string str = to_string(ai[i]);
			glRasterPos2d(cps[i][0] + 30., cps[i][1] + 30.);
			for (char& c : str)
			{
				glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, c);
			}

		}
	}
	else
	{

		for (int i = 0; i < Num - 2; ++i)		//开式ek-curve
		{
			string str = to_string(ai[i]);
			glRasterPos2d(cps[i + 1][0] + 30., cps[i + 1][1] + 30.);
			for (char& c : str)
			{
				glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, c);
			}

		}
	}
		

	if (flag0)
	{
		endcurve();
	}
	
	
	int num1 = (int)MtrxCps_o.size();		

	/*for (int i = 0; i < num1; ++i)
	{
		int n = MtrxCps_o[i].size();
		for (int j = 0; j < n; j++)
		{

			cout << MtrxCps_o[i][j][0] << "," << MtrxCps_o[i][j][1] << ' ';

		}
		cout << endl;

	}*/

	int num2 = (int)MtrxCps_c.size();

	/*for (int i = 0; i < num2; ++i)
	{
		int n = MtrxCps_c[i].size();
		for (int j = 0; j < n; j++)
		{

			cout << MtrxCps_c[i][j][0] << "," << MtrxCps_c[i][j][1] << ' ';

		}
		cout << endl;

	}*/

	for (int i = 0; i < num1; ++i)
	{
		glColor3f(1, 1, 0);		//黄色
		glPointSize(7);
		glBegin(GL_POINTS);			//画点

		int n = MtrxCps_o[i].size();

		for (int j = 0; j < n ; j++)
		{
			glVertex2d(MtrxCps_o[i][j][0], MtrxCps_o[i][j][1]);
		}

		glEnd();
	}

	for (int i = 0; i < num2; ++i)
	{
		glColor3f(1, 1, 0);		//黄色
		glPointSize(7);
		glBegin(GL_POINTS);			//画点

		int n = (int)MtrxCps_c[i].size();

		for (int j = 0; j < n; j++)
		{
			glVertex2d(MtrxCps_c[i][j][0], MtrxCps_c[i][j][1]);
		}

		glEnd();
	}


	//for (int i = 0; i < num1; ++i)
	//{
	//	glColor3f(0, 1, 0);		//绿色
	//	glLineWidth(3);
	//	glBegin(GL_LINE_STRIP);			//画线

	//	int n = MtrxCps_o[i].size();

	//	for (int j = 0; j < n; j++)
	//	{
	//		glVertex2d(MtrxCps_o[i][j][0], MtrxCps_o[i][j][1]);
	//	}
	//	//glVertex2d(MtrxCps[i][0][0], MtrxCps[i][0][1]);
	//	glEnd();
	//}

	//for (int i = 0; i < num2; ++i)
	//{
	//	glColor3f(0, 1, 0);		//绿色
	//	glLineWidth(3);
	//	glBegin(GL_LINE_LOOP);			//画线

	//	int n = MtrxCps_c[i].size();

	//	for (int j = 0; j < n; j++)
	//	{
	//		glVertex2d(MtrxCps_c[i][j][0], MtrxCps_c[i][j][1]);
	//	}
	//	//glVertex2d(MtrxCps[i][0][0], MtrxCps[i][0][1]);
	//	glEnd();
	//}

	for (int i = 0; i < num1; ++i)			//画多条开式ek-curve
	{
		opened_ekcurve(MtrxCps_o[i],Mtrxai_o[i]);
		int n = (int)MtrxCps_o[i].size();
		for (int j = 0; j < n-2; j++)
		{
			drawQuBzr(ci0[j], ci1[j], ci2[j],Mtrxai_o[i][j]);
		}

	}

	for (int i = 0; i < num2; ++i)			//画多条闭式ek-curve
	{
		closed_ekcurve(MtrxCps_c[i],Mtrxai_c[i]);
		int n = (int)MtrxCps_c[i].size();
		for (int j = 0; j < n; j++)
		{
			drawQuBzr(ci0[j], ci1[j], ci2[j],Mtrxai_c[i][j]);
		}

	}

	glFlush();
	//glutSwapBuffers();
	
}
void drawQuBzr(EVec2d &ci0, EVec2d &ci1, EVec2d &ci2,double ai)
{
	int n = 20;
	double t;
	double h = 1. / n;
	EVec2d bezPts;
	glColor3f(1, 0, 0);		//红色
	glLineWidth(3);
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i <= n; ++i)
	{
		t = i * h;
		bezPts = pow(1. - t, 3.0) * ci0 + 3.0 * pow(1. - t, 2.0) * t * ((1.0 - ai) * ci0 + ai * ci1)
			+ 3.0 * (1.0 - t) * pow(t, 2.0) * (ai * ci1 + (1.0 - ai) * ci2) + pow(t, 3.0) * ci2;
		glVertex2d(bezPts[0], bezPts[1]);
	}
	glEnd();
}

void drawCrvtrVctrClmn(EVec2d& ci0, EVec2d& ci1, EVec2d& ci2, double ai)
{
	
	int n = 20;			//画曲率直方图
	vecEg2dd k0;
	k0.resize(n+1);

	double t;
	double h = 1. / n;
	double coeff = 3000.;
	double norm;
	double val0;

	EVec2d ct, dct, ddct, e0, e1;
	EVec2d p0, p1, p2, p3;

	p0 = ci0;
	p1 = (1. - ai) * ci0 + ai * ci1;
	p2 = ai * ci1 + (1. - ai) * ci2;
	p3 = ci2;

	//求曲率的向量和模长
	for (int i = 0; i <= n; ++i) {
		t = i * h;
		dct = 3. * pow(1. - t, 2.) * (p1 - p0) + 6. * (1. - t) * t * (p2 - p1) + 3. * pow(t, 2.) * (p3 - p2);
		ddct = 6. * (1. - t) * (p2 - 2. * p1 + p0) + 6. * t * (p3 - 2. * p2 + p1);
		norm = sqrt(pow(dct[0], 2.) + pow(dct[1], 2.));
		val0 = dct[0] * ddct[1] - dct[1] * ddct[0];
		e0 = dct / norm;
		e1[0] = e0[1];			//e0为切线单位向量，e1为与e0垂直的向量
		e1[1] = -e0[0];			//x0*y0+x1*y1=0, x0=y1,y0=-x1
		k0[i] = coeff * e1 * val0 / pow(norm, 3.);		//k0为曲率乘以系数

	}

	//画曲率端点连线
	glColor3d(0.8549, 0.4392, 0.839215);		//紫色
	glLineWidth(2);
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i <= n; ++i)
	{
		t = i * h;
		ct = pow(1. - t, 3.0) * p0 + 3.0 * pow(1. - t, 2.0) * t * p1
			+ 3.0 * (1.0 - t) * pow(t, 2.0) * p2 + pow(t, 3.0) * p3;
		glVertex2d(ct[0]+k0[i][0], ct[1]+k0[i][1]);
	}
	glEnd();

	//画曲率长度线
	glColor3d(0.8549, 0.4392, 0.839215);		//紫色
	glLineWidth(2);
	glBegin(GL_LINES);
	for (int i = 0; i <= n; ++i)
	{
		t = i * h;
		ct = pow(1. - t, 3.0) * p0 + 3.0 * pow(1. - t, 2.0) * t * p1
			+ 3.0 * (1.0 - t) * pow(t, 2.0) * p2 + pow(t, 3.0) * p3;
		glVertex2d(ct[0], ct[1]);
		glVertex2d(ct[0] + k0[i][0], ct[1] + k0[i][1]);
	}
	glEnd();

}



void myKeyboard(unsigned char key, int x, int y)
{
	if (key == 'e' || key == 'E')
	{
		flag0 = true;
		cout << "e" << endl;
		
		
	}
	glutPostRedisplay();
}


void opened_ekcurve(vecEg2dd& cpsi, vector<double>& ai)
{
	int PntNum = (int)cpsi.size();
	int CntNum = PntNum - 2;
	lambdai.clear();
	lambdai.resize(CntNum, 0.5);
	ti.clear();
	ti.resize(CntNum, 0.5);
	ci0.clear();
	ci0.resize(CntNum);
	ci1.clear();
	ci1.resize(CntNum);
	ci2.clear();
	ci2.resize(CntNum);

	/*ai.clear();
	ai.resize(CntNum, 2. / 3.);*/

	ci0[0] = cpsi[0];
	ci2[CntNum - 1] = cpsi[PntNum - 1];
	for (int i = 0; i < CntNum; ++i)
	{
		ci1[i] = cpsi[i + 1];
	}

	opencalculate_ci02(CntNum);

	int counter = 0;
	const int ITER_NUM = 400;
	
	do{
		/*std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "iter : " << iter << std::endl;*/
		vecEg2dd prevci1;
		prevci1.resize(CntNum);
		prevci1 = ci1;

		opencalculate_ci1(CntNum, cpsi, ai);

		opencalculate_ci02(CntNum);
		opencalculate_ti(CntNum, cpsi, ai);

		opencalculate_lambdai(CntNum, ai);
		++counter;
		if (CheckItr(prevci1, CntNum))
		{
			cout << "converge" << endl;
			break;
		}


	} while (counter < ITER_NUM && !openCheckCvrg(CntNum, cpsi, ai));

}

void closed_ekcurve(vecEg2dd& cpsi, vector<double>& ai)
{
	int numi = (int)cpsi.size();
	lambdai.clear();
	lambdai.resize(numi, 0.5);
	ti.clear();
	ti.resize(numi, 0.5);
	ci0.clear();
	ci0.resize(numi);
	ci1.clear();
	ci1.resize(numi);
	ci2.clear();
	ci2.resize(numi);
	
	/*ai.clear();
	ai.resize(numi, 2. / 3.);*/
	

	ci1 = cpsi;
	calculate_ci02(numi);

	int counter = 0;
	const int ITER_NUM = 400;
	do {
		vecEg2dd prevci1;
		prevci1.resize(numi);
		prevci1 = ci1;

		calculate_ci1(numi, cpsi,ai);
		calculate_ci02(numi);
		calculate_ti(numi, cpsi,ai);
		calculate_lambdai(numi, ai);
		++counter;
		if (CheckItr(prevci1, numi))
		{
			cout << "converge" << endl;
			break;
		}
	} while (counter < ITER_NUM && !closedCheckCvrg(numi, cpsi, ai));
	
}


bool CheckItr(const vecEg2dd &prevci1,int n)
{
	double s = 0;
	for (int i = 0; i < n; ++i)
	{
		s += (ci1[i] - prevci1[i]).norm();
	}
	/*cout << "s=" << s << endl;*/
	return s < 1.e-4;
}

bool openCheckCvrg(int n, vecEg2dd& cpsi, vector<double>& ai)
{
	EVec2d cti;

	double dis0;
	double dis1 = 0;
	for (int i = 0; i < n; ++i)
	{
		cti = pow(1 - ti[i], 3.0) * ci0[i] + 3.0 * pow(1 - ti[i], 2.0) * ti[i] * ((1.0 - ai[i]) * ci0[i] + ai[i] * ci1[i])
			+ 3.0 * (1.0 - ti[i]) * pow(ti[i], 2.0) * (ai[i] * ci1[i] + (1.0 - ai[i]) * ci2[i]) + pow(ti[i], 3.0) * ci2[i];
		
		dis0 = (cti - cpsi[i+1]).norm();
		dis1 += dis0;
	}

	if (dis1 > 1.0e-8)
		return false;

	return true;
}

bool closedCheckCvrg(int n, vecEg2dd& cpsi, vector<double>& ai)
{
	EVec2d cti;
	
	double dis0;
	double dis1 = 0;
	
	for (int i = 0; i < n; ++i)
	{
		cti = pow(1 - ti[i], 3.0) * ci0[i] + 3.0 * pow(1 - ti[i], 2.0) * ti[i] * ((1.0 - ai[i]) * ci0[i] + ai[i] * ci1[i])
			+ 3.0 * (1.0 - ti[i]) * pow(ti[i], 2.0) * (ai[i] * ci1[i] + (1.0 - ai[i]) * ci2[i]) + pow(ti[i], 3.0) * ci2[i];
		
		dis0 = (cti - cpsi[i]).norm();
		dis1 += dis0;
	}

	/*cout << "cti:" << endl;		
	for (EVec2d &ct : cti)
	{
		cout << ct.transpose() << endl;
	}
	cout << "dis0:" << endl;
	for (double &dis : dis0)
	{
		cout << dis << endl;
	}*/
	
	if (dis1 > 1.0e-8)
		return false;
	
	
	return true;
}
void calculate_ti(int n, vecEg2dd& cpsi, vector<double>& ai)
{
	for (int i = 0; i < n; ++i)
	{
		ti[i] = Newton(i, cpsi[i],ai[i]);
		if (ti[i] < 0.0)
		{
			ti[i] = 0.0;
		}
		else if (ti[i] > 1.0)
		{
			ti[i] = 1.0;
		}
	}
	
}

void opencalculate_ti(int n, vecEg2dd& cpsi, vector<double>& ai)
{
	for (int i = 0; i < n; ++i)
	{
		ti[i] = Newton(i, cpsi[i+1],ai[i]);
		if (ti[i] < 0.0)
		{
			ti[i] = 0.0;
		}
		else if (ti[i] > 1.0)
		{
			ti[i] = 1.0;
		}
	}
}

double Newton(int i, EVec2d& cpsii, double aii)
{
	double EPS = 1.e-12;
	double x0, x1;
	x0 = 0;
	x1 = 0.5;
	int k = 0;
	int MaxIter = 100;
	while ((fabs(x1 - x0) >= EPS) && (k<MaxIter))
	{
		x0 = x1;
		x1 = x0 - dCrvtr(i, x0, cpsii, aii) / ddCrvtr(i, x0, cpsii, aii);
		if (k > MaxIter)
		{
			cerr << "Nexton failed \n";
			return 0.;
		}
		k++;
	}
	return x1;
}

//double bisection(int i,EVec2d& cpsii, double aii)
//{
//	double a;					//二分法求曲率的导数为零时ti的值，即曲率极值处的ti值
//	double max = 1.0;
//	double min = 0.0;
//	double val0, val1;
//	val0 = dCrvtr(i, min, cpsii,aii);
//	val1 = dCrvtr(i, max, cpsii,aii);
//
//	if (val0 * val1 >= 0.0)
//	{
//		cerr << "bisection error: f(a)*f(b)>=0 \n";
//		return 0;
//	}
//	double EPS = 1.e-12;
//	while ((max - min) >= EPS)
//	{
//		a = (max + min) / 2.0;
//		if (dCrvtr(i, a, cpsii, aii) == 0.0) { break; }
//		else{
//			if (val0 * dCrvtr(i, a, cpsii, aii) < 0.0) {
//				max = a;
//			}
//			else {
//				min = a;
//			}
//
//		}
//	}
//	return a;
//
//}

double dCrvtr(int i, double ti, EVec2d& cpsii, double aii)
{
	double value;									//使用Maxima求出曲率的导数关于ti的表达式
	double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;

	t0 = -3.0 * pow(ci0[i][1], 2.0) * aii + 6.0 * cpsii[1] * ci0[i][1] * aii - 3.0 * pow(cpsii[1], 2.0) * aii
		- 3.0 * pow(ci0[i][0], 2.0) * aii + 6.0 * cpsii[0] * ci0[i][0] * aii - 3.0 * pow(cpsii[0], 2.0) * aii
		+ 4.0 * pow(ci0[i][1], 2.0) - 8.0 * cpsii[1] * ci0[i][1] + 4.0 * pow(cpsii[1], 2.0) + 4.0 * pow(ci0[i][0], 2.0)
		- 8.0 * cpsii[0] * ci0[i][0] + 4.0 * pow(cpsii[0] , 2.0);

	t1 = (18.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0) - 18.0 * cpsii[1] * ci2[i][1] * pow(aii, 2.0)
		- 18.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		+ 18.0 * cpsii[1] * ci0[i][1] * pow(aii, 2.0)
		+ 18.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		- 18.0 * cpsii[0] * ci2[i][0] * pow(aii, 2.0)
		- 18.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0)
		+ 18.0 * cpsii[0] * ci0[i][0] * pow(aii, 2.0)
		- 36.0 * ci0[i][1] * ci2[i][1] * aii
		+ 36.0 * cpsii[1] * ci2[i][1] * aii
		+ 48.0 * pow(ci0[i][1], 2.0) * aii
		- 60.0 * cpsii[1] * ci0[i][1] * aii
		+ 12.0 * pow(cpsii[1], 2.0) * aii
		- 36.0 * ci0[i][0] * ci2[i][0] * aii
		+ 36.0 * cpsii[0] * ci2[i][0] * aii
		+ 48.0 * pow(ci0[i][0], 2.0) * aii
		- 60.0 * cpsii[0] * ci0[i][0] * aii
		+ 12.0 * pow(cpsii[0], 2.0) * aii + 18.0 * ci0[i][1] * ci2[i][1]
		- 18.0 * cpsii[1] * ci2[i][1] - 30.0 * pow(ci0[i][1], 2.0)
		+ 42.0 * cpsii[1] * ci0[i][1] - 12.0 * pow(cpsii[1], 2.0)
		+ 18.0 * ci0[i][0] * ci2[i][0] - 18.0 * cpsii[0] * ci2[i][0]
		- 30.0 * pow(ci0[i][0], 2.0) + 42.0 * cpsii[0] * ci0[i][0]
		- 12.0 * pow(cpsii[0], 2.0));

	t2 = -144.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0) + 144.0 * cpsii[1] * ci2[i][1] * pow(aii, 2.0)
		+ 144.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		- 144.0 * cpsii[1] * ci0[i][1] * pow(aii, 2.0)
		- 144.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		+ 144.0 * cpsii[0] * ci2[i][0] * pow(aii, 2.0)
		+ 144.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0)
		- 144.0 * cpsii[0] * ci0[i][0] * pow(aii, 2.0)
		+ 258.0 * ci0[i][1] * ci2[i][1] * aii
		- 258.0 * cpsii[1] * ci2[i][1] * aii
		- 276.0 * pow(ci0[i][1], 2) * aii
		+ 294.0 * cpsii[1] * ci0[i][1] * aii
		- 18.0 * pow(cpsii[1], 2.0) * aii
		+ 258.0 * ci0[i][0] * ci2[i][0] * aii
		- 258.0 * cpsii[0] * ci2[i][0] * aii
		- 276.0 * pow(ci0[i][0], 2) * aii
		+ 294.0 * cpsii[0] * ci0[i][0] * aii
		- 18.0 * pow(cpsii[0], 2.0) * aii
		- 114.0 * ci0[i][1] * ci2[i][1]
		+ 114.0 * cpsii[1] * ci2[i][1] + 126.0 * pow(ci0[i][1], 2.0)
		- 138.0 * cpsii[1] * ci0[i][1] + 12.0 * pow(cpsii[1], 2.0)
		- 114.0 * ci0[i][0] * ci2[i][0]
		+ 114.0 * cpsii[0] * ci2[i][0] + 126.0 * pow(ci0[i][0], 2.0)
		- 138.0 * cpsii[0] * ci0[i][0] + 12.0 * pow(cpsii[0], 2.0);

	t3 = 54.0 * pow(ci2[i][1], 2.0) * pow(aii, 3.0) - 108.0 * ci0[i][1] * ci2[i][1] * pow(aii, 3.0)
		+ 54.0 * pow(ci0[i][1], 2.0) * pow(aii, 3.0) + 54.0 * pow(ci2[i][0], 2.0) * pow(aii, 3.0)
		- 108.0 * ci0[i][0] * ci2[i][0] * pow(aii, 3.0)
		+ 54.0 * pow(ci0[i][0], 2.0) * pow(aii, 3.0) - 162.0 * pow(ci2[i][1], 2.0) * pow(aii, 2.0)
		+ 774.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0)
		- 450.0 * cpsii[1] * ci2[i][1] * pow(aii, 2.0)
		- 612.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		+ 450.0 * cpsii[1] * ci0[i][1] * pow(aii, 2.0)
		- 162.0 * pow(ci2[i][0], 2.0) * pow(aii, 2.0)
		+ 774.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		- 450.0 * cpsii[0] * ci2[i][0] * pow(aii, 2.0)
		- 612.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0)
		+ 450.0 * cpsii[0] * ci0[i][0] * pow(aii, 2.0) + 162.0 * pow(ci2[i][1], 2.0) * aii
		- 1032.0 * ci0[i][1] * ci2[i][1] * aii
		+ 708.0 * cpsii[1] * ci2[i][1] * aii + 882.0 * pow(ci0[i][1], 2.0) * aii
		- 732.0 * cpsii[1] * ci0[i][1] * aii + 12.0 * pow(cpsii[1], 2.0) * aii
		+ 162.0 * pow(ci2[i][0], 2.0) * aii - 1032.0 * ci0[i][0] * ci2[i][0] * aii
		+ 708.0 * cpsii[0] * ci2[i][0] * aii + 882.0 * pow(ci0[i][0], 2.0) * aii
		- 732.0 * cpsii[0] * ci0[i][0] * aii + 12.0 * pow(cpsii[0], 2.0) * aii
		- 54.0 * pow(ci2[i][1], 2.0) + 380.0 * ci0[i][1] * ci2[i][1]
		- 272.0 * cpsii[1] * ci2[i][1] - 334.0 * pow(ci0[i][1], 2.0)
		+ 288.0 * cpsii[1] * ci0[i][1] - 8.0 * pow(cpsii[1], 2.0) - 54.0 * pow(ci2[i][0], 2.0)
		+ 380.0 * ci0[i][0] * ci2[i][0] - 272.0 * cpsii[0] * ci2[i][0]
		- 334.0 * pow(ci0[i][0], 2.0) + 288.0 * cpsii[0] * ci0[i][0] - 8.0 * pow(cpsii[0], 2.0);

	t4 = (-405.0 * pow(ci2[i][1], 2.0) * pow(aii, 3.0)) + 810.0 * ci0[i][1] * ci2[i][1] * pow(aii, 3.0)
		- 405.0 * pow(ci0[i][1], 2.0) * pow(aii, 3.0) - 405.0 * pow(ci2[i][0], 2.0) * pow(aii, 3.0)
		+ 810.0 * ci0[i][0] * ci2[i][0] * pow(aii, 3.0)
		- 405.0 * pow(ci0[i][0], 2.0) * pow(aii, 3.0) + 1080.0 * pow(ci2[i][1], 2.0) * pow(aii, 2.0)
		- 2880.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0)
		+ 720.0 * cpsii[1] * ci2[i][1] * pow(aii, 2.0)
		+ 1800.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		- 720.0 * cpsii[1] * ci0[i][1] * pow(aii, 2.0)
		+ 1080.0 * pow(ci2[i][0], 2.0) * pow(aii, 2.0)
		- 2880.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		+ 720.0 * cpsii[0] * ci2[i][0] * pow(aii, 2.0)
		+ 1800.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0)
		- 720.0 * cpsii[0] * ci0[i][0] * pow(aii, 2.0)
		- 945.0 * pow(ci2[i][1], 2.0) * aii
		+ 2910.0 * ci0[i][1] * ci2[i][1] * aii
		- 1020.0 * cpsii[1] * ci2[i][1] * aii
		- 1965.0 * pow(ci0[i][1], 2.0) * aii
		+ 1020.0 * cpsii[1] * ci0[i][1] * aii
		- 945.0 * pow(ci2[i][0], 2.0) * aii
		+ 2910.0 * ci0[i][0] * ci2[i][0] * aii
		- 1020.0 * cpsii[0] * ci2[i][0] * aii
		- 1965.0 * pow(ci0[i][0], 2.0) * aii
		+ 1020.0 * cpsii[0] * ci0[i][0] * aii + 270.0 * pow(ci2[i][1], 2.0)
		- 900.0 * ci0[i][1] * ci2[i][1] + 360.0 * cpsii[1] * ci2[i][1]
		+ 630.0 * pow(ci0[i][1], 2.0) - 360.0 * cpsii[1] * ci0[i][1]
		+ 270.0 * pow(ci2[i][0], 2.0) - 900.0 * ci0[i][0] * ci2[i][0]
		+ 360.0 * cpsii[0] * ci2[i][0] + 630.0 * pow(ci0[i][0], 2.0)
		- 360.0 * cpsii[0] * ci0[i][0];

	t5 = 1296.0 * pow(ci2[i][1], 2.0) * pow(aii, 3.0) - 2592.0 * ci0[i][1] * ci2[i][1] * pow(aii, 3.0)
		+ 1296.0 * pow(ci0[i][1], 2.0) * pow(aii, 3.0) + 1296.0 * pow(ci2[i][0], 2.0) * pow(aii, 3.0)
		- 2592.0 * ci0[i][0] * ci2[i][0] * pow(aii, 3.0)
		+ 1296.0 * pow(ci0[i][0], 2.0) * pow(aii, 3.0) - 3114.0 * pow(ci2[i][1], 2.0) * pow(aii, 2.0)
		+ 6822.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0)
		- 594.0 * cpsii[1] * ci2[i][1] * pow(aii, 2.0)
		- 3708.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		+ 594.0 * cpsii[1] * ci0[i][1] * pow(aii, 2.0)
		- 3114.0 * pow(ci2[i][0], 2.0) * pow(aii, 2.0)
		+ 6822.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		- 594.0 * cpsii[0] * ci2[i][0] * pow(aii, 2.0)
		- 3708.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0)
		+ 594.0 * cpsii[0] * ci0[i][0] * pow(aii, 2.0)
		+ 2454.0 * pow(ci2[i][1], 2.0) * aii
		- 5700.0 * ci0[i][1] * ci2[i][1] * aii
		+ 792.0 * cpsii[1] * ci2[i][1] * aii
		+ 3246.0 * pow(ci0[i][1], 2.0) * aii
		- 792.0 * cpsii[1] * ci0[i][1] * aii
		+ 2454.0 * pow(ci2[i][0], 2.0) * aii
		- 5700.0 * ci0[i][0] * ci2[i][0] * aii
		+ 792.0 * cpsii[0] * ci2[i][0] * aii
		+ 3246.0 * pow(ci0[i][0], 2.0) * aii
		- 792.0 * cpsii[0] * ci0[i][0] * aii - 636 * pow(ci2[i][1], 2.0)
		+ 1536.0 * ci0[i][1] * ci2[i][1] - 264 * cpsii[1] * ci2[i][1]
		- 900.0 * pow(ci0[i][1], 2.0) + 264 * cpsii[1] * ci0[i][1]
		- 636.0 * pow(ci2[i][0], 2.0) + 1536 * ci0[i][0] * ci2[i][0]
		- 264.0 * cpsii[0] * ci2[i][0] - 900 * pow(ci0[i][0], 2.0)
		+ 264.0 * cpsii[0] * ci0[i][0];

	t6 = (-2268.0 * pow(ci2[i][1], 2.0) * pow(aii, 3.0)) + 4536.0 * ci0[i][1] * ci2[i][1] * pow(aii, 3.0)
		- 2268.0 * pow(ci0[i][1], 2.0) * pow(aii, 3.0)
		- 2268.0 * pow(ci2[i][0], 2.0) * pow(aii, 3.0)
		+ 4536.0 * ci0[i][0] * ci2[i][0] * pow(aii, 3.0)
		- 2268.0 * pow(ci0[i][0], 2.0) * pow(aii, 3.0)
		+ 5004.0 * pow(ci2[i][1], 2.0) * pow(aii, 2.0)
		- 10206.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0)
		+ 198.0 * cpsii[1] * ci2[i][1] * pow(aii, 2.0)
		+ 5202.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		- 198.0 * cpsii[1] * ci0[i][1] * pow(aii, 2.0)
		+ 5004.0 * pow(ci2[i][0], 2.0) * pow(aii, 2.0)
		- 10206.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		+ 198.0 * cpsii[0] * ci2[i][0] * pow(aii, 2.0)
		+ 5202.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0)
		- 198.0 * cpsii[0] * ci0[i][0] * pow(aii, 2.0)
		- 3648.0 * pow(ci2[i][1], 2.0) * aii
		+ 7560.0 * ci0[i][1] * ci2[i][1] * aii
		- 264.0 * cpsii[1] * ci2[i][1] * aii
		- 3912.0 * pow(ci0[i][1], 2.0) * aii
		+ 264.0 * cpsii[1] * ci0[i][1] * aii
		- 3648.0 * pow(ci2[i][0], 2.0) * aii
		+ 7560.0 * ci0[i][0] * ci2[i][0] * aii
		- 264.0 * cpsii[0] * ci2[i][0] * aii
		- 3912.0 * pow(ci0[i][0], 2.0) * aii
		+ 264.0 * cpsii[0] * ci0[i][0] * aii + 880.0 * pow(ci2[i][1], 2.0)
		- 1848.0 * ci0[i][1] * ci2[i][1] + 88.0 * cpsii[1] * ci2[i][1]
		+ 968.0 * pow(ci0[i][1], 2.0) - 88.0 * cpsii[1] * ci0[i][1]
		+ 880.0 * pow(ci2[i][0], 2.0) - 1848.0 * ci0[i][0] * ci2[i][0]
		+ 88.0 * cpsii[0] * ci2[i][0] + 968.0 * pow(ci0[i][0], 2.0)
		- 88.0 * cpsii[0] * ci0[i][0];

	t7 = 2268.0 * pow(ci2[i][1], 2.0) * pow(aii, 3.0) - 4536.0 * ci0[i][1] * ci2[i][1] * pow(aii, 3.0)
		+ 2268.0 * pow(ci0[i][1], 2.0) * pow(aii, 3.0) + 2268.0 * pow(ci2[i][0], 2.0) * pow(aii, 3.0)
		- 4536.0 * ci0[i][0] * ci2[i][0] * pow(aii, 3.0)
		+ 2268.0 * pow(ci0[i][0], 2.0) * pow(aii, 3.0) - 4698.0 * pow(ci2[i][1], 2.0) * pow(aii, 2.0)
		+ 9396.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0)
		- 4698.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0) - 4698.0 * pow(ci2[i][0], 2.0) * pow(aii, 2.0)
		+ 9396.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		- 4698.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0) + 3240.0 * pow(ci2[i][1], 2.0) * aii
		- 6480.0 * ci0[i][1] * ci2[i][1] * aii
		+ 3240.0 * pow(ci0[i][1], 2.0) * aii + 3240.0 * pow(ci2[i][0], 2.0) * aii
		- 6480.0 * ci0[i][0] * ci2[i][0] * aii
		+ 3240.0 * pow(ci0[i][0], 2.0) * aii - 744.0 * pow(ci2[i][1], 2.0)
		+ 1488.0 * ci0[i][1] * ci2[i][1] - 744.0 * pow(ci0[i][1], 2.0)
		- 744.0 * pow(ci2[i][0], 2.0) + 1488.0 * ci0[i][0] * ci2[i][0]
		- 744.0 * pow(ci0[i][0], 2.0);

	t8 = (-1215.0 * pow(ci2[i][1], 2.0) * pow(aii, 3.0)) + 2430.0 * ci0[i][1] * ci2[i][1] * pow(aii, 3.0)
		- 1215.0 * pow(ci0[i][1], 2.0) * pow(aii, 3.0)
		- 1215.0 * pow(ci2[i][0], 2.0) * pow(aii, 3.0)
		+ 2430.0 * ci0[i][0] * ci2[i][0] * pow(aii, 3.0)
		- 1215.0 * pow(ci0[i][0], 2.0) * pow(aii, 3.0)
		+ 2430.0 * pow(ci2[i][1], 2.0) * pow(aii, 2.0)
		- 4860.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0)
		+ 2430.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		+ 2430.0 * pow(ci2[i][0], 2.0) * pow(aii, 2.0)
		- 4860.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		+ 2430.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0) - 1620.0 * pow(ci2[i][1], 2.0) * aii
		+ 3240.0 * ci0[i][1] * ci2[i][1] * aii
		- 1620.0 * pow(ci0[i][1], 2.0) * aii - 1620.0 * pow(ci2[i][0], 2.0) * aii
		+ 3240.0 * ci0[i][0] * ci2[i][0] * aii
		- 1620.0 * pow(ci0[i][0], 2.0) * aii + 360.0 * pow(ci2[i][1], 2.0)
		- 720.0 * ci0[i][1] * ci2[i][1] + 360.0 * pow(ci0[i][1], 2.0)
		+ 360.0 * pow(ci2[i][0], 2.0) - 720.0 * ci0[i][0] * ci2[i][0]
		+ 360.0 * pow(ci0[i][0], 2.0);

	t9 = 270.0 * pow(ci2[i][1], 2.0) * pow(aii, 3.0) - 540.0 * ci0[i][1] * ci2[i][1] * pow(aii, 3.0)
		+ 270.0 * pow(ci0[i][1], 2.0) * pow(aii, 3.0)
		+ 270.0 * pow(ci2[i][0], 2.0) * pow(aii, 3.0)
		- 540.0 * ci0[i][0] * ci2[i][0] * pow(aii, 3.0)
		+ 270.0 * pow(ci0[i][0], 2.0) * pow(aii, 3.0)
		- 540.0 * pow(ci2[i][1], 2.0) * pow(aii, 2.0)
		+ 1080.0 * ci0[i][1] * ci2[i][1] * pow(aii, 2.0)
		- 540.0 * pow(ci0[i][1], 2.0) * pow(aii, 2.0)
		- 540.0 * pow(ci2[i][0], 2.0) * pow(aii, 2.0)
		+ 1080.0 * ci0[i][0] * ci2[i][0] * pow(aii, 2.0)
		- 540.0 * pow(ci0[i][0], 2.0) * pow(aii, 2.0) + 360.0 * pow(ci2[i][1], 2.0) * aii
		- 720.0 * ci0[i][1] * ci2[i][1] * aii
		+ 360.0 * pow(ci0[i][1], 2.0) * aii + 360.0 * pow(ci2[i][0], 2.0) * aii
		- 720.0 * ci0[i][0] * ci2[i][0] * aii
		+ 360.0 * pow(ci0[i][0], 2.0) * aii - 80.0 * pow(ci2[i][1], 2.0)
		+ 160.0 * ci0[i][1] * ci2[i][1] - 80.0 * pow(ci0[i][1], 2.0)
		- 80.0 * pow(ci2[i][0], 2.0) + 160.0 * ci0[i][0] * ci2[i][0]
		- 80.0 * pow(ci0[i][0], 2.0);

	value = t0 + ti * t1 + pow(ti, 2.0) * t2 + pow(ti, 3.0) * t3 
		+ pow(ti, 4.0) * t4 + pow(ti, 5.0) * t5 + pow(ti, 6.0) * t6 
		+ pow(ti, 7.0) * t7 + pow(ti, 8.0) * t8 + pow(ti, 9.0) * t9;

	return value;

}

double ddCrvtr(int i, double ti, EVec2d& cpsii, double aii)
{
	double value;									//使用Maxima求出曲率的导数的导数关于ti的表达式
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;

	t0 = 18. * ci0[i][1] * ci2[i][1] * pow(aii,2.) - 18. * cpsii[1] * ci2[i][1] * pow(aii,2.)
		- 18. * pow(ci0[i][1],2.) * pow(aii,2.) + 18. * cpsii[1] * ci0[i][1] * pow(aii,2.)
		+ 18. * ci0[i][0] * ci2[i][0] * pow(aii,2.) - 18. * cpsii[0] * ci2[i][0] * pow(aii,2.)
		- 18. * pow(ci0[i][0] , 2.) * pow(aii,2.) + 18. * cpsii[0] * ci0[i][0] * pow(aii,2.)
		- 36. * ci0[i][1] * ci2[i][1] * aii + 36. * cpsii[1] * ci2[i][1] * aii
		+ 48. * pow(ci0[i][1],2.) * aii - 60. * cpsii[1] * ci0[i][1] * aii + 12. * pow(cpsii[1] , 2.) * aii
		- 36. * ci0[i][0] * ci2[i][0] * aii + 36. * cpsii[0] * ci2[i][0] * aii
		+ 48. * pow(ci0[i][0] , 2.) * aii - 60. * cpsii[0] * ci0[i][0] * aii + 12. * pow(cpsii[0] , 2.) * aii
		+ 18. * ci0[i][1] * ci2[i][1] - 18. * cpsii[1] * ci2[i][1] - 30. * pow(ci0[i][1],2.)
		+ 42. * cpsii[1] * ci0[i][1] - 12. * pow(cpsii[1] , 2.) + 18. * ci0[i][0] * ci2[i][0]
		- 18. * cpsii[0] * ci2[i][0] - 30. * pow(ci0[i][0] , 2.) + 42. * cpsii[0] * ci0[i][0] - 12. * pow(cpsii[0] , 2.);
	t1 = 2. * ((-144. * ci0[i][1] * ci2[i][1] * pow(aii,2.)) + 144. * cpsii[1] * ci2[i][1] * pow(aii,2.)
		+ 144. * pow(ci0[i][1],2.) * pow(aii,2.)
		- 144. * cpsii[1] * ci0[i][1] * pow(aii,2.)
		- 144. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		+ 144. * cpsii[0] * ci2[i][0] * pow(aii,2.)
		+ 144. * pow(ci0[i][0] , 2.) * pow(aii,2.)
		- 144. * cpsii[0] * ci0[i][0] * pow(aii,2.)
		+ 258. * ci0[i][1] * ci2[i][1] * aii
		- 258. * cpsii[1] * ci2[i][1] * aii
		- 276. * pow(ci0[i][1],2.) * aii
		+ 294. * cpsii[1] * ci0[i][1] * aii
		- 18. * pow(cpsii[1] , 2.) * aii
		+ 258. * ci0[i][0] * ci2[i][0] * aii
		- 258. * cpsii[0] * ci2[i][0] * aii
		- 276. * pow(ci0[i][0] , 2.) * aii
		+ 294. * cpsii[0] * ci0[i][0] * aii
		- 18. * pow(cpsii[0] , 2.) * aii
		- 114. * ci0[i][1] * ci2[i][1]
		+ 114. * cpsii[1] * ci2[i][1] + 126. * pow(ci0[i][1],2.)
		- 138. * cpsii[1] * ci0[i][1] + 12. * pow(cpsii[1] , 2.)
		- 114. * ci0[i][0] * ci2[i][0]
		+ 114. * cpsii[0] * ci2[i][0] + 126. * pow(ci0[i][0] , 2.)
		- 138. * cpsii[0] * ci0[i][0] + 12. * pow(cpsii[0] , 2.));
	t2 = 3. * (54. * pow(ci2[i][1] , 2.) * pow(aii , 3.) - 108. * ci0[i][1] * ci2[i][1] * pow(aii , 3.)
		+ 54. * pow(ci0[i][1],2.) * pow(aii , 3.) + 54. * pow(ci2[i][0] , 2.) * pow(aii , 3.)
		- 108. * ci0[i][0] * ci2[i][0] * pow(aii , 3.)
		+ 54. * pow(ci0[i][0] , 2.) * pow(aii , 3.) - 162. * pow(ci2[i][1] , 2.) * pow(aii,2.)
		+ 774. * ci0[i][1] * ci2[i][1] * pow(aii,2.)
		- 450. * cpsii[1] * ci2[i][1] * pow(aii,2.)
		- 612. * pow(ci0[i][1],2.) * pow(aii,2.)
		+ 450. * cpsii[1] * ci0[i][1] * pow(aii,2.)
		- 162. * pow(ci2[i][0] , 2.) * pow(aii,2.)
		+ 774. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		- 450. * cpsii[0] * ci2[i][0] * pow(aii,2.)
		- 612. * pow(ci0[i][0] , 2.) * pow(aii,2.)
		+ 450. * cpsii[0] * ci0[i][0] * pow(aii,2.)
		+ 162. * pow(ci2[i][1] , 2.) * aii
		- 1032. * ci0[i][1] * ci2[i][1] * aii
		+ 708. * cpsii[1] * ci2[i][1] * aii + 882. * pow(ci0[i][1],2.) * aii
		- 732. * cpsii[1] * ci0[i][1] * aii + 12. * pow(cpsii[1] , 2.) * aii
		+ 162. * pow(ci2[i][0] , 2.) * aii
		- 1032. * ci0[i][0] * ci2[i][0] * aii
		+ 708. * cpsii[0] * ci2[i][0] * aii + 882. * pow(ci0[i][0] , 2.) * aii
		- 732. * cpsii[0] * ci0[i][0] * aii + 12. * pow(cpsii[0] , 2.) * aii
		- 54. * pow(ci2[i][1] , 2.) + 380. * ci0[i][1] * ci2[i][1]
		- 272. * cpsii[1] * ci2[i][1] - 334. * pow(ci0[i][1],2.)
		+ 288. * cpsii[1] * ci0[i][1] - 8. * pow(cpsii[1] , 2.) - 54. * pow(ci2[i][0] , 2.)
		+ 380. * ci0[i][0] * ci2[i][0] - 272. * cpsii[0] * ci2[i][0]
		- 334. * pow(ci0[i][0] , 2.) + 288. * cpsii[0] * ci0[i][0]
		- 8. * pow(cpsii[0] , 2.));
	t3 = 4. * ((-405 * pow(ci2[i][1] , 2.) * pow(aii , 3.)) + 810. * ci0[i][1] * ci2[i][1] * pow(aii , 3.)
		- 405. * pow(ci0[i][1],2.) * pow(aii , 3.)
		- 405. * pow(ci2[i][0] , 2.) * pow(aii , 3.)
		+ 810. * ci0[i][0] * ci2[i][0] * pow(aii , 3.)
		- 405. * pow(ci0[i][0] , 2.) * pow(aii , 3.)
		+ 1080. * pow(ci2[i][1] , 2.) * pow(aii,2.)
		- 2880. * ci0[i][1] * ci2[i][1] * pow(aii,2.)
		+ 720. * cpsii[1] * ci2[i][1] * pow(aii,2.)
		+ 1800. * pow(ci0[i][1],2.) * pow(aii,2.)
		- 720. * cpsii[1] * ci0[i][1] * pow(aii,2.)
		+ 1080. * pow(ci2[i][0] , 2.) * pow(aii,2.)
		- 2880. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		+ 720. * cpsii[0] * ci2[i][0] * pow(aii,2.)
		+ 1800. * pow(ci0[i][0] , 2.) * pow(aii,2.)
		- 720. * cpsii[0] * ci0[i][0] * pow(aii,2.)
		- 945. * pow(ci2[i][1] , 2.) * aii
		+ 2910. * ci0[i][1] * ci2[i][1] * aii
		- 1020. * cpsii[1] * ci2[i][1] * aii
		- 1965. * pow(ci0[i][1],2.) * aii
		+ 1020. * cpsii[1] * ci0[i][1] * aii
		- 945. * pow(ci2[i][0] , 2.) * aii
		+ 2910. * ci0[i][0] * ci2[i][0] * aii
		- 1020. * cpsii[0] * ci2[i][0] * aii
		- 1965. * pow(ci0[i][0] , 2.) * aii
		+ 1020. * cpsii[0] * ci0[i][0] * aii + 270. * pow(ci2[i][1] , 2.)
		- 900. * ci0[i][1] * ci2[i][1] + 360. * cpsii[1] * ci2[i][1]
		+ 630. * pow(ci0[i][1],2.) - 360. * cpsii[1] * ci0[i][1]
		+ 270. * pow(ci2[i][0] , 2.) - 900. * ci0[i][0] * ci2[i][0]
		+ 360. * cpsii[0] * ci2[i][0] + 630. * pow(ci0[i][0] , 2.)
		- 360. * cpsii[0] * ci0[i][0]);
	t4 = 5. * (1296. * pow(ci2[i][1] , 2.) * pow(aii , 3.) - 2592 * ci0[i][1] * ci2[i][1] * pow(aii , 3.)
		+ 1296. * pow(ci0[i][1],2.) * pow(aii , 3.)
		+ 1296. * pow(ci2[i][0] , 2.) * pow(aii , 3.)
		- 2592. * ci0[i][0] * ci2[i][0] * pow(aii , 3.)
		+ 1296. * pow(ci0[i][0] , 2.) * pow(aii , 3.)
		- 3114. * pow(ci2[i][1] , 2.) * pow(aii,2.)
		+ 6822. * ci0[i][1] * ci2[i][1] * pow(aii,2.)
		- 594. * cpsii[1] * ci2[i][1] * pow(aii,2.)
		- 3708. * pow(ci0[i][1],2.) * pow(aii,2.)
		+ 594. * cpsii[1] * ci0[i][1] * pow(aii,2.)
		- 3114. * pow(ci2[i][0] , 2.) * pow(aii,2.)
		+ 6822. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		- 594. * cpsii[0] * ci2[i][0] * pow(aii,2.)
		- 3708. * pow(ci0[i][0] , 2.) * pow(aii,2.)
		+ 594. * cpsii[0] * ci0[i][0] * pow(aii,2.)
		+ 2454. * pow(ci2[i][1] , 2.) * aii
		- 5700. * ci0[i][1] * ci2[i][1] * aii
		+ 792. * cpsii[1] * ci2[i][1] * aii
		+ 3246. * pow(ci0[i][1],2.) * aii
		- 792. * cpsii[1] * ci0[i][1] * aii
		+ 2454. * pow(ci2[i][0] , 2.) * aii
		- 5700. * ci0[i][0] * ci2[i][0] * aii
		+ 792. * cpsii[0] * ci2[i][0] * aii
		+ 3246. * pow(ci0[i][0] , 2.) * aii
		- 792. * cpsii[0] * ci0[i][0] * aii - 636. * pow(ci2[i][1] , 2.)
		+ 1536. * ci0[i][1] * ci2[i][1] - 264. * cpsii[1] * ci2[i][1]
		- 900. * pow(ci0[i][1],2.) + 264. * cpsii[1] * ci0[i][1]
		- 636. * pow(ci2[i][0] , 2.) + 1536. * ci0[i][0] * ci2[i][0]
		- 264. * cpsii[0] * ci2[i][0] - 900. * pow(ci0[i][0] , 2.)
		+ 264. * cpsii[0] * ci0[i][0]);
	t5 = 6. * ((-2268. * pow(ci2[i][1] , 2.) * pow(aii , 3.)) + 4536. * ci0[i][1] * ci2[i][1] * pow(aii , 3.)
		- 2268. * pow(ci0[i][1],2.) * pow(aii , 3.)
		- 2268. * pow(ci2[i][0] , 2.) * pow(aii , 3.)
		+ 4536. * ci0[i][0] * ci2[i][0] * pow(aii , 3.)
		- 2268. * pow(ci0[i][0] , 2.) * pow(aii , 3.)
		+ 5004. * pow(ci2[i][1] , 2.) * pow(aii,2.)
		- 10206. * ci0[i][1] * ci2[i][1] * pow(aii,2.)
		+ 198. * cpsii[1] * ci2[i][1] * pow(aii,2.)
		+ 5202. * pow(ci0[i][1],2.) * pow(aii,2.)
		- 198. * cpsii[1] * ci0[i][1] * pow(aii,2.)
		+ 5004. * pow(ci2[i][0] , 2.) * pow(aii,2.)
		- 10206. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		+ 198. * cpsii[0] * ci2[i][0] * pow(aii,2.)
		+ 5202. * pow(ci0[i][0] , 2.) * pow(aii,2.)
		- 198. * cpsii[0] * ci0[i][0] * pow(aii,2.)
		- 3648. * pow(ci2[i][1] , 2.) * aii
		+ 7560. * ci0[i][1] * ci2[i][1] * aii
		- 264. * cpsii[1] * ci2[i][1] * aii
		- 3912. * pow(ci0[i][1],2.) * aii
		+ 264. * cpsii[1] * ci0[i][1] * aii
		- 3648. * pow(ci2[i][0] , 2.) * aii
		+ 7560. * ci0[i][0] * ci2[i][0] * aii
		- 264. * cpsii[0] * ci2[i][0] * aii
		- 3912. * pow(ci0[i][0] , 2.) * aii
		+ 264. * cpsii[0] * ci0[i][0] * aii + 880. * pow(ci2[i][1] , 2.)
		- 1848. * ci0[i][1] * ci2[i][1]
		+ 88. * cpsii[1] * ci2[i][1] + 968. * pow(ci0[i][1],2.)
		- 88. * cpsii[1] * ci0[i][1] + 880. * pow(ci2[i][0] , 2.)
		- 1848. * ci0[i][0] * ci2[i][0]
		+ 88. * cpsii[0] * ci2[i][0] + 968. * pow(ci0[i][0] , 2.)
		- 88. * cpsii[0] * ci0[i][0]);
	t6 = 7. * (2268. * pow(ci2[i][1] , 2.) * pow(aii , 3.) - 4536. * ci0[i][1] * ci2[i][1] * pow(aii , 3.)
		+ 2268. * pow(ci0[i][1],2.) * pow(aii , 3.)
		+ 2268. * pow(ci2[i][0] , 2.) * pow(aii , 3.)
		- 4536. * ci0[i][0] * ci2[i][0] * pow(aii , 3.)
		+ 2268. * pow(ci0[i][0] , 2.) * pow(aii , 3.)
		- 4698. * pow(ci2[i][1] , 2.) * pow(aii,2.)
		+ 9396. * ci0[i][1] * ci2[i][1] * pow(aii,2.)
		- 4698. * pow(ci0[i][1],2.) * pow(aii,2.)
		- 4698. * pow(ci2[i][0] , 2.) * pow(aii,2.)
		+ 9396. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		- 4698. * pow(ci0[i][0] , 2.) * pow(aii,2.) + 3240. * pow(ci2[i][1] , 2.) * aii
		- 6480. * ci0[i][1] * ci2[i][1] * aii
		+ 3240. * pow(ci0[i][1],2.) * aii + 3240. * pow(ci2[i][0] , 2.) * aii
		- 6480. * ci0[i][0] * ci2[i][0] * aii
		+ 3240. * pow(ci0[i][0] , 2.) * aii - 744. * pow(ci2[i][1] , 2.)
		+ 1488. * ci0[i][1] * ci2[i][1] - 744. * pow(ci0[i][1],2.)
		- 744. * pow(ci2[i][0] , 2.) + 1488. * ci0[i][0] * ci2[i][0]
		- 744. * pow(ci0[i][0] , 2.));
	t7 = 8. * ((-1215. * pow(ci2[i][1] , 2.) * pow(aii , 3.)) + 2430. * ci0[i][1] * ci2[i][1] * pow(aii , 3.)
		- 1215. * pow(ci0[i][1],2.) * pow(aii , 3.)
		- 1215. * pow(ci2[i][0] , 2.) * pow(aii , 3.)
		+ 2430. * ci0[i][0] * ci2[i][0] * pow(aii , 3.)
		- 1215. * pow(ci0[i][0] , 2.) * pow(aii , 3.)
		+ 2430. * pow(ci2[i][1] , 2.) * pow(aii,2.)
		- 4860. * ci0[i][1] * ci2[i][1] * pow(aii,2.)
		+ 2430. * pow(ci0[i][1],2.) * pow(aii,2.)
		+ 2430. * pow(ci2[i][0] , 2.) * pow(aii,2.)
		- 4860. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		+ 2430. * pow(ci0[i][0] , 2.) * pow(aii,2.)
		- 1620. * pow(ci2[i][1] , 2.) * aii
		+ 3240. * ci0[i][1] * ci2[i][1] * aii
		- 1620. * pow(ci0[i][1],2.) * aii - 1620. * pow(ci2[i][0] , 2.) * aii
		+ 3240. * ci0[i][0] * ci2[i][0] * aii
		- 1620. * pow(ci0[i][0] , 2.) * aii + 360. * pow(ci2[i][1] , 2.)
		- 720. * ci0[i][1] * ci2[i][1] + 360. * pow(ci0[i][1],2.)
		+ 360. * pow(ci2[i][0] , 2.) - 720. * ci0[i][0] * ci2[i][0]
		+ 360. * pow(ci0[i][0] , 2.));
	t8 = 9. * (270. * pow(ci2[i][1] , 2.) * pow(aii , 3.) - 540 * ci0[i][1] * ci2[i][1] * pow(aii , 3.)
		+ 270. * pow(ci0[i][1],2.) * pow(aii , 3.)
		+ 270. * pow(ci2[i][0] , 2.) * pow(aii , 3.)
		- 540. * ci0[i][0] * ci2[i][0] * pow(aii , 3.)
		+ 270. * pow(ci0[i][0] , 2.) * pow(aii , 3.)
		- 540. * pow(ci2[i][1] , 2.) * pow(aii,2.)
		+ 1080. * ci0[i][1] * ci2[i][1] * pow(aii,2.)
		- 540. * pow(ci0[i][1],2.) * pow(aii,2.)
		- 540. * pow(ci2[i][0] , 2.) * pow(aii,2.)
		+ 1080. * ci0[i][0] * ci2[i][0] * pow(aii,2.)
		- 540. * pow(ci0[i][0] , 2.) * pow(aii,2.)
		+ 360. * pow(ci2[i][1] , 2.) * aii
		- 720. * ci0[i][1] * ci2[i][1] * aii
		+ 360. * pow(ci0[i][1],2.) * aii + 360. * pow(ci2[i][0] , 2.) * aii
		- 720. * ci0[i][0] * ci2[i][0] * aii
		+ 360. * pow(ci0[i][0] , 2.) * aii - 80. * pow(ci2[i][1] , 2.)
		+ 160. * ci0[i][1] * ci2[i][1] - 80. * pow(ci0[i][1],2.)
		- 80. * pow(ci2[i][0] , 2.) + 160. * ci0[i][0] * ci2[i][0]
		- 80. * pow(ci0[i][0] , 2.));

	value = t0 + ti * t1 + pow(ti, 2.0) * t2 + pow(ti, 3.0) * t3
		+ pow(ti, 4.0) * t4 + pow(ti, 5.0) * t5 + pow(ti, 6.0) * t6
		+ pow(ti, 7.0) * t7 + pow(ti, 8.0) * t8;

	return value;
}


void calculate_lambdai(int n, vector<double>& ai)
{

	for (int i = 0; i < n; ++i)
	{
		
		int nextI = (i + 1) % n;
		double A = (1.0 - ai[i]) * TrgArea(ci0[i], ci1[i], ci1[nextI]);
		double B = (1.0 - ai[nextI]) * TrgArea(ci1[i], ci1[nextI], ci2[nextI]);
		lambdai[i] = sqrt(A) / (sqrt(A)+ (ai[i] / ai[nextI]) * sqrt(B));
		if (lambdai[i] < 0.0)
		{
			lambdai[i] = 0.0;
		}
		else if (lambdai[i] > 1.0)
		{
			lambdai[i] = 1.0;
		}

	}

}

void opencalculate_lambdai(int n, vector<double>& ai)
{
	for (int i = 0; i < n-1; i++)
	{
		double A = (1. - ai[i]) * TrgArea(ci0[i], ci1[i], ci1[i + 1]);
		double B = (1. - ai[i + 1]) * TrgArea(ci1[i], ci1[i + 1], ci2[i + 1]);
		lambdai[i] = sqrt(A) / (sqrt(A) + (ai[i] / ai[i + 1]) * sqrt(B));
		if (lambdai[i] < 0.0)
		{
			lambdai[i] = 0.0;
		}
		else if (lambdai[i] > 1.0)
		{
			lambdai[i] = 1.0;
		}
	}
}

double TrgArea(const EVec2d&p0,const EVec2d&p1,const EVec2d&p2)
{
	EVec2d v1 = p1 - p0;
	EVec2d v2 = p2 - p0;
	double Area = 0.5*fabs(v1[0] * v2[1] - v1[1] * v2[0]);
	return Area;
}

void opencalculate_ci1(int n, vecEg2dd& cpsi, vector<double>& ai)
{
	MatrixXd A(n, n);
	VectorXd b1(n), b2(n);
	A.setZero();
	if (n == 1)
	{
		ci1[0] = (cpsi[1] - pow(1. - ti[0], 3.) * cpsi[0] - 3. * (1. - ti[0]) * ti[0] * (1. - ai[0]) * ((1. - ti[0]) * cpsi[0] + ti[0] * cpsi[2]) - pow(ti[0], 3.) * cpsi[2])
			/ (3. * ai[0] * (1. - ti[0]) * ti[0]);
		return;
	}
	
	for (int i = 0; i < n; ++i)
	{
		if (i == 0)
		{
			A(0, 0) = 3. * pow(1. - ti[0], 2.) * ti[0] * ai[0] 
				+ 3. * (1. - ti[0]) * pow(ti[0], 2.) * (ai[0] + (1. - ai[0]) * (1. - lambdai[0])) 
				+ pow(ti[0], 3.) * (1. - lambdai[0]);
			A(0, 1) = 3. * (1. - ti[0]) * pow(ti[0], 2.) * (1. - ai[0]) * lambdai[0] 
				+ pow(ti[0], 3.) * lambdai[0];
			b1[0] = cpsi[1][0] - pow(1.0 - ti[0], 3.) * cpsi[0][0] - 3. * pow(1. - ti[0], 2.) * ti[0] * (1. - ai[0]) * cpsi[0][0];
			b2[0] = cpsi[1][1] - pow(1.0 - ti[0], 3.) * cpsi[0][1] - 3. * pow(1. - ti[0], 2.) * ti[0] * (1. - ai[0]) * cpsi[0][1];
		}
		else if (i == n - 1)
		{
			A(n - 1, n - 2) = pow(1. - ti[n - 1], 3.) * (1. - lambdai[n - 2]) 
				+ 3. * pow(1. - ti[n - 1], 2.) * ti[n - 1] * (1. - ai[n - 1]) * (1. - lambdai[n - 2]);
			A(n - 1, n - 1) = pow(1. - ti[n - 1], 3.) * lambdai[n - 2] 
				+ 3. * pow(1. - ti[n - 1], 2.) * ti[n - 1] * ((1. - ai[n - 1]) * lambdai[n - 2] + ai[n - 1]) 
				+ 3. * (1. - ti[n - 1]) * pow(ti[n - 1], 2.) * ai[n - 1];
			b1[n - 1] = cpsi[n][0] - (3. * (1. - ti[n - 1]) * pow(ti[n - 1], 2.) * (1. - ai[n - 1]) + pow(ti[n - 1], 3.)) * cpsi[n + 1][0];
			b2[n - 1] = cpsi[n][1] - (3. * (1. - ti[n - 1]) * pow(ti[n - 1], 2.) * (1. - ai[n - 1]) + pow(ti[n - 1], 3.)) * cpsi[n + 1][1];
		}
		else
		{
			A(i, i - 1) = pow(1. - ti[i], 3.) * (1. - lambdai[i - 1]) + 3. * pow(1. - ti[i], 2.) * ti[i] * (1. - ai[i]) * (1. - lambdai[i - 1]);
			A(i, i) = pow(1. - ti[i], 3.) * lambdai[i - 1] + 3. * pow(1. - ti[i], 2.) * ti[i] * ((1. - ai[i]) * lambdai[i - 1] + ai[i])
				+ 3. * (1. - ti[i]) * pow(ti[i], 2.) * ((1. - ai[i]) * (1. - lambdai[i]) + ai[i])
				+ pow(ti[i], 3.) * (1. - lambdai[i]);
			A(i, i + 1) = 3. * (1. - ti[i]) * pow(ti[i], 2.) * (1. - ai[i]) * lambdai[i] + pow(ti[i], 3.) * lambdai[i];
			b1[i] = cpsi[i + 1][0];
			b2[i] = cpsi[i + 1][1];
		}
	}
	VectorXd x1 = A.lu().solve(b1);
	VectorXd x2 = A.lu().solve(b2);

	for (int i = 0; i < n; ++i)
	{
		ci1[i][0] = x1[i];
		ci1[i][1] = x2[i];

	}
}

void calculate_ci1(int n,vecEg2dd &cpsi, vector<double>& ai)
{
	MatrixXd A(n, n);
	VectorXd b1(n), b2(n);
	A.setZero();
	for (int i = 0; i < n; ++i)
	{
		int nextI = (i + 1) % n;
		int prevI = (i - 1+ n) % n;
		A(i, prevI) = pow(1.0 - ti[i], 3.0) * (1.0 - lambdai[prevI]) + 3.0 * pow(1.0 - ti[i], 2.0) * ti[i] * (1.0 - ai[i]) * (1.0 - lambdai[prevI]);
		A(i, i) = pow(1.0 - ti[i], 3.0) * lambdai[prevI] + 3.0 * pow(1.0 - ti[i], 2.0) * ti[i] * ((1.0 - ai[i]) * lambdai[prevI] + ai[i]) + 3.0 * (1.0 - ti[i]) * pow(ti[i], 2.0) * ((1.0 - ai[i]) * (1.0 - lambdai[i]) + ai[i]) + pow(ti[i], 3.0) * (1.0 - lambdai[i]);
		A(i, nextI) = 3.0 * (1.0 - ti[i]) * pow(ti[i], 2.0) * (1.0 - ai[i]) * lambdai[i] + pow(ti[i], 3.0) * lambdai[i];
		b1[i] = cpsi[i][0];
		b2[i] = cpsi[i][1];
	}
	/*cout << "Here is the invertible matrix A:" << endl << A << endl;*/
	
	VectorXd x1 = A.lu().solve(b1);
	VectorXd x2 = A.lu().solve(b2);
	
	for (int i = 0; i < n; ++i)
	{
		ci1[i][0] = x1[i];
		ci1[i][1] = x2[i];
	}
	/*cout << "Here is the matrix ci1:" << endl;
	for (EVec2d &point : ci1)
	{
		cout << point.transpose() << endl;
	}*/

}

void opencalculate_ci02(int n)
{
	for (int i = 0; i < (n-1); ++i)
	{
		ci2[i] = ci0[i + 1] = (1 - lambdai[i]) * ci1[i] + lambdai[i] * ci1[i + 1];
	}
}

void calculate_ci02(int n)
{
	for (int i = 0; i < n; ++i)
	{
		int nextI = (i + 1) % n;
		ci2[i] = ci0[nextI] = (1 - lambdai[i])*ci1[i] + lambdai[i] * ci1[nextI];
	}
	/*cout << "Here is the matrix ci2:" << endl;
	for (EVec2d &point : ci2)
	{
		cout << point.transpose() << endl;
	}
	cout << "Here is the matrix ci0:" << endl;
	for (EVec2d &point : ci0)
	{
		cout << point.transpose() << endl;
	}*/
}

void changeai(int i)
{
	double a;

	if (closedflag)
	{
		cout << "Please enter the value of a(2/3<=a<1)" << endl;
		cin >> a;
		if (a >= 2. / 3. && a < 1)
		{
			ai[i]=a;
		}
		else {
			ai[i]=defaultai;
		}

	}
	else
	{
		if ( cps.size() > 1)
		{
			cout << "Please enter the value of a(2/3<=a<1)" << endl;
			cin >> a;
			if (a >= 2. / 3. && a < 1)
			{
				ai[i-1]=a;
			}
			else {
				ai[i-1]=defaultai;
			}
		}
	}
}

int getaiIndex(int x, int y)
{
	for (int i = 0; i< Num; i++)
	{
		
		if (fabs(cps[i][0] - (double)x) < 5. || fabs(cps[i][1] - (double)y) < 5.)
			return i;

	}
	return -1;
}

void OnMouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		/*if (Num != 0 && checkpt(x, (winHeight - y)))
		{
			mouseLeftDown = true;
		}
		else
		{
			cps.push_back(EVec2d(x, winHeight - y));
						
		}*/
		
		if (!aiflag)
		{
			cps.push_back(EVec2d(x, winHeight - y));

			if ((Num>3) &&(checkpt(x, winHeight - y)))
			{
				cps[0][0] = x;					//若终点与起点重合，则起始点合为一点，变为闭式ek-curve
				cps[0][1] = winHeight - y;
				cps.pop_back();
				closedflag = true;
				

			}
			
			
		}
		else
		{
			int i = getaiIndex(x, winHeight - y);
			if (i >= 0) {
				changeai(i);
			}
			
			aiflag = false;
		}
		
		
		
	}
	/*else if (state == GLUT_UP)
	{
		mouseLeftDown = false;
	}*/

	glutPostRedisplay();

}

void MouseMove(int x, int y)
{

	cps[Num - 1][0] = x;
	cps[Num - 1][1] = winHeight - y;
	
	glutPostRedisplay();
}

bool checkpt(int x, int y)
{
		
	if (fabs(cps[0][0] - (double)x) < 10. && fabs(cps[0][1] - (double)y) < 10.)
	{
		return true;
	}
		

	return false;
}
//void mouseMotionPT(int x, int y)
//{
//	if (mouseLeftDown)
//	{
//
//		cps[i0][0] = (double)x;
//		cps[i0][1] = (double)(winHeight - y);
//		glutPostRedisplay();
//	}
//}

void menuFunc(int value)
{
	switch (value)
	{
	case 1:
		closedflag = false;
		break;

	case 2:
		aiflag = true;
		
		break;

	case 999:
		exit(1);
		break;
	}
	glutPostRedisplay();
}

void saveandload(int value)
{
	switch (value)
	{
	case 3:
		SvPtsFl();
		break;
	case 4:
		InpPtsFl();
		break;
	
	}
	glutPostRedisplay();
}

void SvPtsFl()
{
	// The filter for txt files
	const char* filePatterns[1] = { "*.txt" };

	// Get the file name
	const char* filename = tinyfd_saveFileDialog("Save Ek curves", "Ekcurve.txt", 1, filePatterns, NULL);

	ofstream outfile(filename);

	if (!outfile.is_open())
		cout << "Open file failed" << endl;
	if (filename)
	{
		//存储封闭ek-curve控制点
		int num1 = (int)Mtrxai_c.size();
		for (int i = 0; i < num1; ++i)
		{
			// Get the number of Bezier curves
			int Pntnum = (int)MtrxCps_c[i].size();

			// Write the number of bezier curves found
			outfile << "arcs <" << Pntnum << ">" << '\t' << "# number of arcs following" << '\n';
			outfile << "arc <" << Pntnum << ">" << '\t' << "# degree = number of control points - 1" << '\n';

			outfile << "closed" << '\t' << "# type of ek-curve" << '\n';

			// Traverse through the control points of the curve
			for (int j = 0; j < Pntnum; j++ )
			{
				
					//		// Write the current control point (only x and y)
				outfile << "<" << MtrxCps_c[i][j].x() << " " << MtrxCps_c[i][j].y() << " " << Mtrxai_c[i][j] << ">" << '\t' << "# control point as two floats separated by blank" << '\n';
			}

			//	// Write the end of the arc
			outfile << "endarc" << '\t' << "# terminates arc description" << '\n';
		}

		//存储开式ek-curve控制点
		int num2 = (int)Mtrxai_o.size();
		for (int i = 0; i < num2; ++i)
		{
			// Get the number of Bezier curves
			int Pntnum = (int)MtrxCps_o[i].size();

			// Write the number of bezier curves found
			outfile << "arcs <" << Pntnum << ">" << '\t' << "# number of arcs following" << '\n';
			outfile << "arc <" << Pntnum << ">" << '\t' << "# degree = number of control points - 1" << '\n';

			outfile << "opened" << '\t' << "# type of ek-curve" << '\n';

			// Traverse through the control points of the curve
			for (int j = 0; j < Pntnum; j++)
			{

				//		// Write the current control point (only x and y)
				outfile << "<" << MtrxCps_o[i][j].x() << " " << MtrxCps_o[i][j].y() << " " << Mtrxai_o[i][j] << ">" << '\t' << "# control point as two floats separated by blank" << '\n';
			}

			//	// Write the end of the arc
			outfile << "endarc" << '\t' << "# terminates arc description" << '\n';
		}
		
	}

	// Close the file
	outfile.close();
	
}

void InpPtsFl()
{
	// The filter for txt files
	const char* filePatterns[1] = { "*.txt" };

	// Display the open file dialog
	const char* filename = tinyfd_openFileDialog("Import a Points File", "", 1, filePatterns, NULL, 0);

	// Open the file as an input stream file 
	ifstream infile(filename);
	// The string where individual lines will be stored
	string line;

	if (!infile.is_open())
		cout << "Open file failed" << endl;

	if (filename)
	{
		// Read the file line by line
		while (getline(infile, line))
		{
			// Split the line with the comment char
			vector<string> tokens = Utils::split(line, '#');

			// Discard the comments, split by blanks
			tokens = Utils::split(tokens.at(0), ' ');

			// Take the command of the line and trim it
			string first = tokens.at(0);
			Utils::trim_inplace(first);

			// Compare the command and check which one is
			if (first.compare("arcs") == 0)
			{
				// DO NOTHING!
				// This field is irrelevant since curves are allocated dinamycally

				std::cout << "arcs = " << tokens.at(1) << std::endl;
			}
			else if (first.compare("arc") == 0)
			{
				// DO NOTHING!
				// This field is irrelevant since curve points are allocated dinamycally

				std::cout << "arc = " << tokens.at(1) << std::endl;
			}

			else if (first.compare("closed") == 0)
			{
				closedflag = true;
				cout << "closed ek-curve" << endl;
			}

			else if (first.compare("opened") == 0)
			{
				closedflag = false;
				cout << "opened ek-curve" << endl;
			}

			else if (first.compare("endarc") == 0)
			{
				std::cout << "endarc" << std::endl;
				endcurve();
			}

			else
			{
				// Get the X value of the point (trim blanks and <>)
				std::string sx = tokens.at(0);
				Utils::trim_inplace(sx);
				Utils::trim_inplace(sx, "<>");
				double x = std::stod(sx);

				// Get the Y value of the point (trim blanks and <>)
				std::string sy = tokens.at(1);
				Utils::trim_inplace(sy);
				Utils::trim_inplace(sy, "<>");
				double y = std::stod(sy);

				// Get the a value
				std::string sa = tokens.at(2);
				Utils::trim_inplace(sa);
				Utils::trim_inplace(sa, "<>");
				double a = std::stod(sa);

				std::cout << "double (" << x << ", " << y << "," << a << ")" << std::endl;

				// Generate a new point ans push it into the control points for the new current curve
				cps.push_back(EVec2d(x, y));
				ai.push_back(a);

			}

		}
	}

	infile.close();

}

int main(int argc, char **argv)
{
	int menu,submenu1;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(winwidth, winHeight);
	glutInitWindowPosition(400, 150);
	glutCreateWindow("ek-curve");
	myInit();
	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);
	glutMouseFunc(OnMouse);
	glutMotionFunc(MouseMove);
	glutKeyboardFunc(myKeyboard);


	submenu1 = glutCreateMenu(saveandload);
	glutAddMenuEntry("Save", 3);
	glutAddMenuEntry("Import", 4);
	menu = glutCreateMenu(menuFunc);
	glutAddMenuEntry("Draw Ek-curve", 1);
	glutAddMenuEntry("Change a", 2);
	glutAddSubMenu("Save & Import", submenu1);
	glutAddMenuEntry("Exit", 999);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	//glutMotionFunc(mouseMotionPT);
	glutMainLoop();
	return 0;
}