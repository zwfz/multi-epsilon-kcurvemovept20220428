#include<gl/glut.h>
#include<vector>
#include<Eigen/Dense>
#include<iostream>
#include<math.h>

//#include <glew.h>
//#include <freeglut.h>
//#include <freeglut_ext.h>
/*使用三次方程求根法求ti*/
/*19/1/26 计算连接处的切线角度，给LAC环赋角度初值*/

using namespace std;
using namespace Eigen;

typedef Eigen::Vector2d EVec2d;
typedef vector<EVec2d, Eigen::aligned_allocator<EVec2d>> vecEg2dd;

typedef vector<vecEg2dd, Eigen::aligned_allocator<vecEg2dd>> VvecEg2dd;			//二维点列的矩阵容器

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

GLsizei winwidth = 800;
GLsizei winHeight = 800;
bool mouseLeftDown;
bool mouseRightDown;
int Num,i0;
bool flag0 = false;
bool closedflag = true;

void kcurve(vecEg2dd& cpsi);
void openkcurve(vecEg2dd& cpsi);

void closed_ekcurve(vecEg2dd& cpsi);

void calculate_ci02(int n);
void opencalculate_ci02(int n);
void calculate_ci1(int n, vecEg2dd& cpsi);
void opencalculate_ci1(int n, vecEg2dd& cpsi);
void calculate_ti(int n, vecEg2dd& cpsi);
void opencalculate_ti(int n, vecEg2dd& cpsi);
void calculate_lambdai(int n);
void opencalculate_lambdai(int n);
double TrgArea(const EVec2d&p0, const EVec2d&p1, const EVec2d&p2);
bool CheckCvrg(int n, vecEg2dd& cpsi);
bool openCheckCvrg(int n, vecEg2dd& cpsi);
bool CheckItr(const vecEg2dd &prevci1,int n);
void drawQuBzr(EVec2d &ci0, EVec2d &ci1, EVec2d &ci2);
//bool checkpt(int x, int y);
//void mouseMotionPT(int x, int y);
void OnMouse(int button, int state, int x, int y);
void myKeyboard(unsigned char key, int x, int y);

void menuFunc(int value);
//void multikcrv();
double getRealSolutionOfCubicFunc(const double a, const double b, const double c);

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

void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	Num = cps.size();

	for (EVec2d& point : cps)
	{
		cout << "cps:" << point.transpose() << endl;
	}

	glColor3f(1, 1, 0);		//黄色
	glPointSize(7);
	glBegin(GL_POINTS);

	for (int i = 0; i < Num; ++i)
	{
		glVertex2d(cps[i][0], cps[i][1]);
	}

	glEnd();

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
			kcurve(cps);
			for (int j = 0; j < Num; ++j)
			{
				drawQuBzr(ci0[j], ci1[j], ci2[j]);
			}
		}
		if (!closedflag)
		{
			openkcurve(cps);
			for (int j = 0; j < Num-2; ++j)
			{
				drawQuBzr(ci0[j], ci1[j], ci2[j]);
			}
		}
		
	}
		

	if (flag0)
	{
		if (closedflag)
		{
			MtrxCps_c.push_back({ cps });			//闭式K-curve点矩阵
		}
		if (!closedflag)
		{
			MtrxCps_o.push_back({ cps });			//开式K-curve点矩阵
		}
		
		flag0 = false;
		cps.clear();
	}
	

	
	int num1 = MtrxCps_o.size();		

	for (int i = 0; i < num1; ++i)
	{
		int n = MtrxCps_o[i].size();
		for (int j = 0; j < n; j++)
		{

			cout << MtrxCps_o[i][j][0] << "," << MtrxCps_o[i][j][1] << ' ';

		}
		cout << endl;

	}

	int num2 = MtrxCps_c.size();

	for (int i = 0; i < num2; ++i)
	{
		int n = MtrxCps_c[i].size();
		for (int j = 0; j < n; j++)
		{

			cout << MtrxCps_c[i][j][0] << "," << MtrxCps_c[i][j][1] << ' ';

		}
		cout << endl;

	}

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

		int n = MtrxCps_c[i].size();

		for (int j = 0; j < n; j++)
		{
			glVertex2d(MtrxCps_c[i][j][0], MtrxCps_c[i][j][1]);
		}

		glEnd();
	}


	for (int i = 0; i < num1; ++i)
	{
		glColor3f(1, 0, 1);		//粉红色
		glLineWidth(3);
		glBegin(GL_LINE_STRIP);			//画线

		int n = MtrxCps_o[i].size();

		for (int j = 0; j < n; j++)
		{
			glVertex2d(MtrxCps_o[i][j][0], MtrxCps_o[i][j][1]);
		}
		//glVertex2d(MtrxCps[i][0][0], MtrxCps[i][0][1]);
		glEnd();
	}

	for (int i = 0; i < num2; ++i)
	{
		glColor3f(1, 0, 1);		//粉红色
		glLineWidth(3);
		glBegin(GL_LINE_LOOP);			//画线

		int n = MtrxCps_c[i].size();

		for (int j = 0; j < n; j++)
		{
			glVertex2d(MtrxCps_c[i][j][0], MtrxCps_c[i][j][1]);
		}
		//glVertex2d(MtrxCps[i][0][0], MtrxCps[i][0][1]);
		glEnd();
	}

	for (int i = 0; i < num1; ++i)			//画多条开式k-curve
	{
		openkcurve(MtrxCps_o[i]);
		int n = MtrxCps_o[i].size();
		for (int j = 0; j < n-2; j++)
		{
			drawQuBzr(ci0[j], ci1[j], ci2[j]);
		}

	}

	for (int i = 0; i < num2; ++i)			//画多条闭式k-curve
	{
		kcurve(MtrxCps_c[i]);
		int n = MtrxCps_c[i].size();
		for (int j = 0; j < n; j++)
		{
			drawQuBzr(ci0[j], ci1[j], ci2[j]);
		}

	}

	glFlush();
	//glutSwapBuffers();
	
}
void drawQuBzr(EVec2d &ci0, EVec2d &ci1, EVec2d &ci2)
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
		bezPts = (1. - t)*(1. - t)*ci0 + 2.*(1. - t)*t*ci1 + t * t*ci2;
		glVertex2d(bezPts[0], bezPts[1]);
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

//void multikcrv() {
//	cout << "multikcrv" << endl;
//	MtrxCps.push_back({ cps });
//	int num1 = MtrxCps.size();
//	for (int i = 0; i < num1; ++i)
//	{
//		Num = cps.size();
//		cout << "Num:" << Num << endl;
//		cout << "i:" << i << endl;
//		//if(Num>3)
//		kcurve(MtrxCps[i]);
//	}
//	cps.clear();
//}

void openkcurve(vecEg2dd& cpsi)
{
	int PntNum = cpsi.size();
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
	cti.clear();
	cti.resize(CntNum);
	dcti.clear();
	dcti.resize(CntNum);
	dcps.clear();
	dcps.resize(CntNum);
	theta.clear();
	theta.resize(CntNum);

	ci0[0] = cpsi[0];
	ci2[CntNum - 1] = cpsi[PntNum - 1];
	for (int i = 0; i < CntNum; ++i)
	{
		ci1[i] = cpsi[i + 1];
	}
	
	opencalculate_ci02(CntNum);
	

	const int ITER_NUM = 400;
	for (int iter = 0; iter < ITER_NUM; ++iter)
	{
		/*std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "iter : " << iter << std::endl;*/
		vecEg2dd prevci1;
		prevci1.resize(CntNum);
		prevci1 = ci1;

		opencalculate_ci1(CntNum, cpsi);

		opencalculate_ci02(CntNum);
		opencalculate_ti(CntNum, cpsi);

		opencalculate_lambdai(CntNum);

		if (openCheckCvrg(CntNum, cpsi) || CheckItr(prevci1, CntNum))
		{
			cout << "converge" << endl;
			break;
		}

	}

}


void closed_ekcurve(vecEg2dd& cpsi)
{
	int numi = cpsi.size();
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
	cti.clear();
	cti.resize(numi);
	dcti.clear();
	dcti.resize(numi);
	dcps.clear();
	dcps.resize(numi);
	theta.clear();
	theta.resize(numi);

	ci1 = cpsi;
	calculate_ci02(numi);

	int counter = 0;
	const int ITER_NUM = 400;
	do {
		vecEg2dd prevci1;
		prevci1.resize(numi);
		prevci1 = ci1;

		calculate_ci1(numi, cpsi);
		calculate_ci02(numi);
		calculate_ti(numi, cpsi);
	}
}

void kcurve(vecEg2dd &cpsi)
{
	int numi = cpsi.size();
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
	cti.clear();
	cti.resize(numi);
	dcti.clear();
	dcti.resize(numi);
	dcps.clear();
	dcps.resize(numi);
	theta.clear();
	theta.resize(numi);

	ci1 = cpsi;
	calculate_ci02(numi);

	const int ITER_NUM = 400;
	/*cout << "ITER_NUM="<<ITER_NUM << endl;*/
	for (int iter = 0; iter < ITER_NUM; ++iter)
	{
		/*std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "iter : " << iter << std::endl;*/
		vecEg2dd prevci1;
		prevci1.resize(numi);
		prevci1 = ci1;

		calculate_ci1(numi,cpsi);
		
		calculate_ci02(numi);
		calculate_ti(numi, cpsi);
		
		calculate_lambdai(numi);
		
		if (CheckCvrg(numi, cpsi) || CheckItr(prevci1, numi))
		{
			cout<<"converge"<<endl;
			break;
		}
			
	}
	

	/*glColor3d(0, 1, 1);
	glLineWidth(3);
	glBegin(GL_LINE_LOOP);

	for (int i = 0; i < Num; ++i)
	{

		glVertex2d(ci0[i][0], ci0[i][1]);
		glVertex2d(ci1[i][0], ci1[i][1]);
	}
	glEnd();*/

	/*for (int i = 0; i < numi; ++i)
		drawQuBzr(ci0[i], ci1[i], ci2[i]);*/

	/*求控制点连线的向量*/
	for (int i = 0; i < numi; ++i)
	{
		dcps[i] = cpsi[(i + 1)%numi] - cpsi[i];
	}
	cout << "cps:" << endl;
	for (EVec2d &ct : cpsi)
	{
		cout << ct.transpose() << endl;
	}
	cout << "dcps:" << endl;
	for (EVec2d &ct : dcps)
	{
		cout << ct.transpose() << endl;
	}

	//glColor3f(0, 1, 0);			//绿色
	//glPointSize(5);
	//glBegin(GL_POINTS);
	//for (int i = 0; i < numi; ++i)
	//{
	//	
	//	cti[i] = (1 - ti[i])*(1 - ti[i])*ci0[i] + 2 * (1 - ti[i])*ti[i] * ci1[i] + ti[i] * ti[i] * ci2[i];
	//	dcti[i] = 2 * ((1 - ti[i])*(ci1[i] - ci0[i]) + ti[i] * (ci2[i] - ci1[i]));	/*求控制点处切线的向量*/
	//	glVertex2d(cti[i][0], cti[i][1]);
	//}
	//glEnd();
	//
	///*求控制点切线与k曲线在控制点处切线的夹角*/
	//for (int i = 0; i < numi; ++i)
	//{
	//	double s0,s1;
	//	s0 = dcti[i][0] * dcps[i][0] + dcti[i][1] * dcps[i][1];
	//	s1 = sqrt(dcps[i][0] * dcps[i][0] + dcps[i][1] * dcps[i][1])
	//		 *sqrt(dcti[i][0] * dcti[i][0] + dcti[i][1] * dcti[i][1]);
	//	theta[i] = acos(s0 / s1)*180/PI;
	//	
	//}
	//cout << "theta:" << endl;
	//for (double &ct : theta)
	//{
	//	cout << ct << endl;
	//}

	//glColor3f(0, 1, 1);			//青色
	//glLineWidth(3);
	//for (int i = 0; i < numi; ++i)
	//{
	//	EVec2d dctf,dctb;
	//	dctf[0] = cti[i][0]+100.0;
	//	dctf[1] = cti[i][1] + 100.0 * dcti[i][1] / dcti[i][0];
	//	dctb[0] = cti[i][0] - 100.0;
	//	dctb[1] = cti[i][1] - 100.0 * dcti[i][1] / dcti[i][0];
	//	glBegin(GL_LINE_STRIP);
	//	glVertex2d(dctb[0], dctb[1]);
	//	glVertex2d(dctf[0], dctf[1]);
	//	glEnd();
	//}
	

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

bool openCheckCvrg(int n, vecEg2dd& cpsi)
{
	vecEg2dd cti;
	cti.resize(n);
	vector<double> dis0;
	dis0.resize(n);
	for (int i = 0; i < n; ++i)
	{
		cti[i][0] = (1 - ti[i]) * (1 - ti[i]) * ci0[i][0] + 2 * (1 - ti[i]) * ti[i] * ci1[i][0] + ti[i] * ti[i] * ci2[i][0];
		cti[i][1] = (1 - ti[i]) * (1 - ti[i]) * ci0[i][1] + 2 * (1 - ti[i]) * ti[i] * ci1[i][1] + ti[i] * ti[i] * ci2[i][1];
		dis0[i] = (cti[i] - cpsi[i+1]).norm();
	}

	for (int i = 0; i < n; ++i)
	{
		if (dis0[i] > 1.0e-4)
			return false;
	}

	return true;
}

bool CheckCvrg(int n, vecEg2dd& cpsi)
{
	vecEg2dd cti;
	cti.resize(n);
	vector<double> dis0;
	dis0.resize(n);
	for (int i = 0; i < n; ++i)
	{
		cti[i][0] = (1 - ti[i])*(1 - ti[i])*ci0[i][0] + 2 * (1 - ti[i])*ti[i] * ci1[i][0] + ti[i] * ti[i] * ci2[i][0];
		cti[i][1] = (1 - ti[i])*(1 - ti[i])*ci0[i][1] + 2 * (1 - ti[i])*ti[i] * ci1[i][1] + ti[i] * ti[i] * ci2[i][1];
		dis0[i] = (cti[i] - cpsi[i]).norm();
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
	for (int i = 0; i < n; ++i)
	{
		if (dis0[i] > 1.0e-4)
			return false;
	}
	
	return true;
}
void calculate_ti(int n, vecEg2dd& cpsi)
{
	for (int i = 0; i < n; ++i)
	{
		EVec2d ci2_ci0 = ci2[i] - ci0[i];
		EVec2d ci0_pi = ci0[i] - cpsi[i];
		double a = ci2_ci0.squaredNorm();
		double b = 3 * ci2_ci0.dot(ci0_pi);
		double c = (3 * ci0[i] - 2 * cpsi[i] - ci2[i]).dot(ci0_pi);
		double d = -ci0_pi.squaredNorm();
		ti[i] = (a == 0) ? 0.5 : getRealSolutionOfCubicFunc(b / a, c / a, d / a);
	}
	
}

void opencalculate_ti(int n, vecEg2dd& cpsi)
{
	for (int i = 0; i < n; ++i)
	{
		EVec2d ci2_ci0 = ci2[i] - ci0[i];
		EVec2d ci0_piplus1 = ci0[i] - cpsi[i + 1];
		double a = ci2_ci0.squaredNorm();
		double b = 3 * ci2_ci0.dot(ci0_piplus1);
		double c = (3 * ci0[i] - 2 * cpsi[i+1] - ci2[i]).dot(ci0_piplus1);
		double d = -ci0_piplus1.squaredNorm();
		ti[i] = (a == 0) ? 0.5 : getRealSolutionOfCubicFunc(b / a, c / a, d / a);
	}
}

double getRealSolutionOfCubicFunc(const double a, const double b, const double c)
{
	const double p = b - a * a / 3;
	const double q = 2 * a*a*a / 27 - a * b / 3 + c;
	const double tmp = q * q / 4 + p * p*p / 27;
	if (tmp < 0) return 0;
	const double mTmp = -q / 2 + sqrt(tmp);
	const double nTmp = -q / 2 - sqrt(tmp);
	const double m = (mTmp < 0) ? -pow(-mTmp, 1 / 3.0) : pow(mTmp, 1 / 3.0);
	const double n = (nTmp < 0) ? -pow(-nTmp, 1 / 3.0) : pow(nTmp, 1 / 3.0);
	return -a / 3 + m + n;
}
void calculate_lambdai(int n)
{
	vector<double> prevLambdai = lambdai;
	for (int i = 0; i < n; ++i)
	{
		
		int nextI = (i + 1) % n;
		double A = sqrt(TrgArea(ci0[i], ci1[i], ci1[nextI]));
		double B = sqrt(TrgArea(ci1[i], ci1[nextI], ci2[nextI]));
		lambdai[i] = A / (A + B);
	}
}

void opencalculate_lambdai(int n)
{
	for (int i = 0; i < n-1; i++)
	{
		double A = sqrt(TrgArea(ci0[i],ci1[i],ci1[i+1]));
		double B = sqrt(TrgArea(ci1[i],ci1[i+1],ci2[i+1]));
		lambdai[i] = A / (A + B);
	}
}

double TrgArea(const EVec2d&p0,const EVec2d&p1,const EVec2d&p2)
{
	EVec2d v1 = p1 - p0;
	EVec2d v2 = p2 - p0;
	double Area = 0.5*fabs(v1[0] * v2[1] - v1[1] * v2[0]);
	return Area;
}

void opencalculate_ci1(int n, vecEg2dd& cpsi)
{
	MatrixXd A(n, n);
	VectorXd b1(n), b2(n);
	A.setZero();
	if (n == 1)
	{
		ci1[0] = (cpsi[1] - pow(1 - ti[0], 2) * cpsi[0] - pow(ti[0], 2) * cpsi[2]) / (2 * ti[0] * (1 - ti[0]));
		return;
	}
	
	for (int i = 0; i < n; ++i)
	{
		if (i == 0)
		{
			A(0, 0) = 2.0 * ti[0] * (1.0 - ti[0]) + ti[0] * ti[0] * (1.0 - lambdai[0]);
			A(0, 1) = ti[0] * ti[0] * lambdai[0];
			b1[0] = cpsi[1][0] - pow(1.0 - ti[0], 2) * cpsi[0][0];
			b2[0] = cpsi[1][1] - pow(1.0 - ti[0], 2) * cpsi[0][1];
		}
		else if (i == n - 1)
		{
			A(n - 1, n - 2) = pow(1.0 - ti[n - 1], 2) * (1.0 - lambdai[n - 2]);
			A(n - 1, n - 1) = pow(1.0 - ti[n - 1], 2) * lambdai[n - 2] + 2.0 * ti[n - 1] * (1.0 - ti[n - 1]);
			b1[n - 1] = cpsi[n][0] - ti[n - 1] * ti[n - 1] * cpsi[n + 1][0];
			b2[n - 1] = cpsi[n][1] - ti[n - 1] * ti[n - 1] * cpsi[n + 1][1];
		}
		else
		{
			A(i, i - 1) = (1.0 - lambdai[i - 1]) * pow(1.0 - ti[i], 2);
			A(i, i) = lambdai[i - 1] * pow(1.0 - ti[i], 2) + (2.0 - (1.0 + lambdai[i]) * ti[i]) * ti[i];
			A(i, i + 1) = lambdai[i] * ti[i] * ti[i];
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

void calculate_ci1(int n,vecEg2dd &cpsi)
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
		cps.push_back(EVec2d(x, winHeight - y));
		
	}
	/*else if (state == GLUT_UP)
	{
		mouseLeftDown = false;
	}*/

	glutPostRedisplay();

}
//bool checkpt(int x, int y)
//{
//	for (i0 = 0; i0 < Num; i0++)
//	{
//		
//		if (fabs(cps[i0][0] - (double)x) < 2. || fabs(cps[i0][1] - (double)y) < 2.)
//			return true;
//
//	}
//	return false;
//}
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
		closedflag = true;
		break;

	case 999:
		exit(1);
		break;
	}
}

int main(int argc, char **argv)
{
	int menu;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(winwidth, winHeight);
	glutInitWindowPosition(400, 150);
	glutCreateWindow("ek-curve");
	myInit();
	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);
	glutMouseFunc(OnMouse);
	glutKeyboardFunc(myKeyboard);

	menu = glutCreateMenu(menuFunc);
	glutAddMenuEntry("Opened Ek-curve", 1);
	glutAddMenuEntry("Closed Ek-curve", 2);
	glutAddMenuEntry("Exit", 999);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	//glutMotionFunc(mouseMotionPT);
	glutMainLoop();
	return 0;
}