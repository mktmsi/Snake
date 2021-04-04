//12-9-1
#include <iostream>
#include <stdio.h>
#include<math.h>
#include <cmath>
#include <stdlib.h>
#include <GLUT/glut.h>
#include <time.h>
#include "Monitor.hpp"
#include "Vector.h"
#include <string.h>
#include <png.h>
#define dt (0.001)	//タイムステップ0.0001
#define MAX_STEP (30000)//最大タイムステップ数
#define PI (3.14159265359)
#define N (30)//質点数
#define trunk(i,j) trunk[N*(j)+(i)]
#define leg(i,j) leg[N*(j)+(i)]
#define r(i,j) r[N*(j)+(i)]
#define v(i,j) v[N*(j)+(i)]

Monitor monitor;


typedef struct{Vector2D r,v;} AA;
AA operator +(AA a,AA b){AA c;	c.r=a.r+b.r;	c.v=a.v+b.v;	return c;}
AA operator -(AA a,AA b){AA c;	c.r=a.r-b.r;	c.v=a.v-b.v ;	return c;}
AA ax(double a,AA b){AA c;	c.r=b.r*a;		c.v=b.v*a;	return c;}
AA wx(double a,AA b){AA c;	c.r=b.r/a;		c.v=b.v/a;	return c;}


//AA trunk[6*N];

/*********************質点********************/
class Mass{
private:
public:
    double m,battery;
    //Vector2D F;
    AA aa;
    //int state_num;
};

Mass trunk[6*(N+1)];
Mass leg[6*(N+1)];

//質点の位置・相互関係
Vector2D vecttor_tt[N],vecttor_tl[N];//質点間の方向ベクトル
Vector2D et_tt[N],et_tl[N];//質点間の単位方向ベクトル
Vector2D en_tt[N],en_tl[N];//質点間の単位法線ベクトル
Vector2D et_tt_hat[N],et_tl_hat[N];//体軸方向の単位ベクトル
Vector2D en_tt_hat[N],en_tl_hat[N];//体軸方向に対して垂直な単位方向ベクトル
double abs_tt[N],abs_tl[N];//質点間の距離
double l_ttDot[N],l_tlDot[N];//質点間距離の変化量
double l_tbar=10.0,l_lbar=2.0;//体節間の自然長
double l_t0=l_tbar,l_l0=l_lbar;//体節間の初期長，脚の書記長
double phi_t[N],phi_tDot[N],phi_tBar[N];//i番目の質点の角度
double phi_l[N],phi_lDot[N],phi_lBar[N];//i番目の質点の角度
//質点にかかる力・トルク
Vector2D F_tBody[N],F_lBody[N];//バネ・ダンパ由来の力ベクトル
Vector2D F_tTorque[N],F_lTorque[N];//トルク由来の力ベクトル
Vector2D F_tFriction[N],F_lFriction[N];//摩擦抵抗の力ベクトル
double ft_s[N],ft_d[N],fl_s[N],fl_d[N];
double tau_t[N],tau_tPas[N],tau_tAct[N],tau_l[N],tau_lPas[N],tau_lAct[N];

//機械定数
double kt=1000.0,kt_pas=5.1,kt_act=20.0;//バネ定数kpas=15.1,kt_act=50.0
double kl=1000.0,kl_pas=5.1,kl_act=20.0;
double ct=100.5,ct_pas=1.11;//ダンパ定数cpas=1.1
double cl=100.5,cl_pas=1.11;//ダンパ定数cpas=1.1
double mt=0.001,mn=1.0;//摩擦定数 mn=0.001
double m=5.0;//質量5.0

//先頭関節のトルクa*sin(wt)
double a=1.0;//a=2.0
double w=20.0;//100.0;

double t;
double L=100.0;
int winid;
int ts;//現在のタイムステップ
int initial=1;
FILE *fp;	//fpというファイル用変数を定義


void capture(int *);

const int save_flag=0;    //0:画像保存しない，1:画像保存する

void init(){
  int i,j;
  double theta;
  char filename[256];
  t,ts=0;
  //等方性摩擦
  //mt=mn;
  sprintf(filename, "data.txt");    //数値データ保存ファイルのファイル名設定
  fp=fopen(filename,"w");        //数値データ保存ファイルを開く

  printf("initial place was initialized\n");
  for(i=0;i<N;i++){
    trunk(i,0).aa.r.set_x(-l_t0*i);
    trunk(i,0).aa.r.set_y(l_l0);//r(i,0).set_x();
    leg(i,0).aa.r.set_x(-l_t0*i);
    leg(i,0).aa.r.set_y(0.0);//r(i,0).set_x();   
  }
  printf("x0_x=%lf,v0_x=%lf\n",trunk[0].aa.r.get_x(),trunk[0].aa.v.get_x());
}


void func(Mass *pt, Mass *pt_out,Mass *pl, Mass *pl_out){
  int i;
  Vector2D temp;
  double kaku;
  
  //質点の位置関係に関する変数を更新
  for(i=1;i<N;i++){
    //体節間の相対関係
    vecttor_tt[i] = pt[i-1].aa.r - pt[i].aa.r;//質点間の相対距離
    abs_tt[i] = vecttor_tt[i].get_abs();
    et_tt[i] = vecttor_tt[i] / abs_tt[i];
    en_tt[i].set_x(et_tt[i].get_y()*sin(-PI/2));//回転行列
    en_tt[i].set_y(et_tt[i].get_x()*sin(PI/2));//回転行列
    l_ttDot[i] = (pt[i-1].aa.v-pt[i].aa.v) * et_tt[i];//体軸方向の相対速度性分

  }

//体軸方向の単位ベクトルを更新
  for(i=1;i<N-1;i++){
    temp = et_tt[i] + et_tt[i+1];
    et_tt_hat[i] = temp/temp.get_abs();
    en_tt_hat[i].set_x(et_tt_hat[i].get_y()*sin(-PI/2));//回転行列
    en_tt_hat[i].set_y(et_tt_hat[i].get_x()*sin(PI/2));//回転行列
  }
  //i=0のときの体軸方向の単位ベクトルを更新
  et_tt_hat[0].set_x(et_tt_hat[1].get_x()*cos(phi_t[1])+et_tt_hat[1].get_y()*sin(-phi_t[1]));//回転行列
  et_tt_hat[0].set_y(et_tt_hat[1].get_x()*sin(phi_t[1])+et_tt_hat[1].get_y()*cos(phi_t[1]));//回転行列
  en_tt_hat[0].set_x(et_tt_hat[0].get_y()*sin(-PI/2));//回転行列
  en_tt_hat[0].set_y(et_tt_hat[0].get_x()*sin(PI/2));//回転行列

  //i=N-1のときの体軸方向の単位ベクトルを更新
  et_tt_hat[N-1].set_x(et_tt_hat[N-2].get_x()*cos(-phi_t[N-2])+et_tt_hat[N-2].get_y()*sin(phi_t[N-2]));//回転行列
  et_tt_hat[N-1].set_y(et_tt_hat[N-2].get_x()*sin(-phi_t[N-2])+et_tt_hat[N-2].get_y()*cos(-phi_t[N-2]));//回転行列
  en_tt_hat[N-1].set_x(et_tt_hat[N-1].get_y()*sin(-PI/2));//回転行列
  en_tt_hat[N-1].set_y(et_tt_hat[N-1].get_x()*sin(PI/2));//回転行列



  //各質点にかかるバネ・ダンパ由来の力ベクトルを更新
  for(i=1;i<N;i++){
    ft_s[i] = kt*pow(abs_tt[i]-l_tbar,3);//spring
    ft_d[i] = ct*l_ttDot[i];//damper
  }
  for(i=0;i<N;i++){
    if(i==0){
      F_tBody[i] = et_tt[i+1]*(-(ft_s[i+1]+ft_d[i+1]));
    }else if(i==N-1){
      F_tBody[i] = et_tt[i]*(ft_s[i]+ft_d[i]);
    }else{
      F_tBody[i] = et_tt[i]*(ft_s[i]+ft_d[i])+et_tt[i+1]*(-(ft_s[i+1]+ft_d[i+1]));
    }
  }

  //質点の角度を更新
  for(i=1;i<N;i++){
    if(i==N-1){
      phi_t[i] = atan2(pt[i-1].aa.r.get_y()-pt[i].aa.r.get_y(),pt[i-1].aa.r.get_x()-pt[i].aa.r.get_x());
      phi_tDot[i] = (pt[i-1].aa.v-pt[i].aa.v)*en_tt[i]/abs_tt[i];
    }else{
      phi_t[i] = atan2(et_tt[i+1]^et_tt[i],et_tt[i+1]*et_tt[i]);
      phi_tDot[i] = ((pt[i-1].aa.v-pt[i].aa.v)*en_tt[i])/abs_tt[i] + ((pt[i+1].aa.v-pt[i].aa.v)*en_tt[i+1])/abs_tt[i+1];
    }
  }

  //目標角度・質点にかかるトルクを更新
  for(i=1;i<N-1;i++){
    if(i==1){
      tau_t[i] = a*sin(w*t);
      tau_t[i] += -kt_pas*phi_t[i] - ct_pas*phi_tDot[i];
    }else{
      phi_tBar[i] = phi_t[i-1];
      tau_tPas[i] = -kt_pas*phi_t[i] - ct_pas*phi_tDot[i];
      tau_tAct[i] = -kt_act*(phi_t[i] - phi_tBar[i]);
      tau_t[i] = tau_tPas[i] + tau_tAct[i];
    }
  }

  //質点にかかるトルク由来の力ベクトルを更新
  for(i=0;i<N;i++){
    if(i==0){
      F_tTorque[0] = en_tt[1]*(tau_t[1]/abs_tt[1]);
    }else if(i==1){
      F_tTorque[1] = en_tt[1]*(-tau_t[1]/abs_tt[1]) + en_tt[2]*((tau_t[2]-tau_t[1])/abs_tt[2]);
    }else if(i==N-2){
      F_tTorque[N-2] = en_tt[N-2]*((tau_t[N-3]-tau_t[N-2])/abs_tt[N-2]) + en_tt[N-1]*(-tau_t[N-2]/abs_tt[N-1]);
    }else if(i==N-1){
      F_tTorque[N-1] = en_tt[N-1]*(tau_t[N-2]/abs_tt[N-1]);
    }else{
      F_tTorque[i] = en_tt[i]*((tau_t[i-1]-tau_t[i])/abs_tt[i]) + en_tt[i+1]*((tau_t[i+1]-tau_t[i])/abs_tt[i+1]);
    }
  }

  //摩擦力を更新
  for(i=0;i<N;i++){
    F_tFriction[i] = et_tt_hat[i]*((-mt)*(pt[i].aa.v*et_tt_hat[i])) + en_tt_hat[i]*((-mn)*(pt[i].aa.v*en_tt_hat[i]));
  }

  //運動方程式
  for(i=0;i<N;i++){
    pt_out[i].aa.v=(F_tBody[i] + F_tTorque[i] + F_tFriction[i])*(dt/m);
    pt_out[i].aa.r = pt[i].aa.v*dt;
  }
}


void runge(){
  int i,j,q,n;
  //i番目の質点の位置に関する微分方程式を解く
  for (q=0; q<100; q++) {
    //printf("runge  x0_x=%lf,v0_x=%lf\n",trunk(0,0).aa.r.get_x(),trunk(0,0).aa.v.get_x());
    //UpdateUnitVectors(&trunk(0,0).aa);
    //printf("runge2  x0_x=%lf,v0_x=%lf\n",trunk(0,0).aa.r.get_x(),trunk(0,0).aa.v.get_x());
    func(&trunk(0,0),&trunk(0,1),&leg(0,0),&leg(0,1));    //trunk(0,1)=k1
    for (j=0; j<N; j++) {
      trunk(j,5).aa=trunk(j,0).aa+wx(2.0,trunk(j,1).aa);
      leg(j,5).aa=leg(j,0).aa+wx(2.0,leg(j,1).aa);
    }
    //UpdateUnitVectors(&trunk(0,5).aa);
    func(&trunk(0,5),&trunk(0,2),&leg(0,5),&leg(0,2));    //trunk(0,2)=k2
    for (j=0; j<N; j++) {
      trunk(j,5).aa=trunk(j,0).aa+wx(2.0,trunk(j,2).aa);
      leg(j,5).aa=leg(j,0).aa+wx(2.0,leg(j,2).aa);
    }
    //UpdateUnitVectors(&trunk(0,5).aa);
    func(&trunk(0,5),&trunk(0,3),&leg(0,5),&leg(0,3));    //trunk(0,3)=k3
    for (j=0; j<N; j++) {
      trunk(j,5).aa=trunk(j,0).aa+trunk(j,3).aa;
      leg(j,5).aa=leg(j,0).aa+leg(j,3).aa;
    }
    //UpdateUnitVectors(&trunk(0,5).aa);
    func(&trunk(0,5),&trunk(0,4),&leg(0,5),&leg(0,4));    //trunk(0,4)=k4
    for (j=0; j<N; j++) {
      trunk(j,0).aa=trunk(j,0).aa+wx(6.0,trunk(j,1).aa+ax(2.0,trunk(j,2).aa)+ax(2.0,trunk(j,3).aa)+trunk(j,4).aa); //x(t+dt)=x(t)+k
      leg(j,0).aa=leg(j,0).aa+wx(6.0,leg(j,1).aa+ax(2.0,leg(j,2).aa)+ax(2.0,leg(j,3).aa)+leg(j,4).aa); //x(t+dt)=x(t)+k
    }
  }

  ts++;
  t=ts*dt;
  //結果の表示
  if(ts%100==0){
    //printf("t=%f  \n",(double)t);
    if(ts>= MAX_STEP){ //MAX_STEPに達したらプログラムを終了
      fclose(fp);
      exit(0);
    }
  }
}



void keyboard(unsigned char key, int x , int y){
  double tmpcin;
  switch(key)
  {
    case 'Q'   :                                break;
    case 'q'   :    exit(0);                            break;
    case '\033':  /* '\033' = ESC */ exit(0);   break;
    case 'h'   : monitor.SetCenter( 1.0 , 0 );  break;
    case 'l'   : monitor.SetCenter( -1.0 , 0 ); break;
    case 'j'   : monitor.SetCenter( 0 , 1.0 );  break;
    case 'k'   : monitor.SetCenter( 0 , -1.0 ); break;
    case 'c'   :                                break;
    case 'g'   :                                break;
    case 'z'   : monitor.SetZoom( 1.1/ 1.0 );   break;
    case 'x'   : monitor.SetZoom( 1.0/ 1.1 );  break;
  }

}

void idle(void)
{
  glutSetWindow(winid);
  //glutKeyboardFunc(keyboard);
  glutPostRedisplay();
}

void display(void)
{
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  int i,j;
  double o=0.0,n;
  //ディスプレイ表示用文字列
  //char str[256],str1[256],str2[256],str3[256],str4[256],str5[256];
  char str[256];
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  //glutKeyboardFunc(keyboard);

  runge();        //runge-kuttaを回す
  //printf("x0_x=%lf,v0_x=%lf\n",trunk[0].aa.r.get_x(),trunk[0].aa.v.get_x());
  /*ィンドウサイズに表示*/
  monitor.SetCenter(trunk(j,0).aa.r.get_x()/L*2-1.5  , 0 );
  //  o=n;
  for(j=0;j<N;j++){
    if(tau_t[j]>0){//トルクが正であればその質点は青色
      monitor.SetAllColor(0.0,0.0,1.0);//青
    }else{
      monitor.SetAllColor(1.0,0.0,.00);//赤
    }
    monitor.DrawCircle(trunk(j,0).aa.r.get_x()/L*2-1.5,trunk(j,0).aa.r.get_y()/L*2,0.02);
  }

  monitor.SetAllColor(0.0,0.0,0.0);
  for(j=0;j<N;j++){
    if(j!=N-1){
      monitor.DrawLine(trunk(j,0).aa.r.get_x()/L*2-1.5,trunk(j,0).aa.r.get_y()/L*2,trunk(j+1,0).aa.r.get_x()/L*2-1.5,trunk(j+1,0).aa.r.get_y()/L*2,0.005);

    }
  }
  for(j=-N/3;j<300;j++){
    monitor.DrawLine(0.3*j,1,0.3*(j),-2,0.005);
  }
  monitor.SetAllColor(0.0,0.0,0.0);

  double x=0.1,y=0.6,dy=0.05;
  sprintf(str,"Isotropic friction");
  monitor.String(x,y,str);
     sprintf(str, "ts=%d", ts);
     y = y - dy;
    monitor.String(x, y, str);
    y = y - dy;
    sprintf(str, "N=%d", N);
    monitor.String(x, y, str);
    y = y - dy;


  glFlush();
  glutSwapBuffers();

  //1000ステップごとに画像を保存
  if(!(ts%100)) {
    //capture(&ts);
  }
}




void mouse(int button, int state, int x, int y){
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) { glutIdleFunc(0);     std::cout << "left: on" << std::endl;  }
    else                    { glutIdleFunc(idle);  std::cout << "left: off" << std::endl; }
    break;

    case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN) { glutIdleFunc(0);    std::cout << "middle: on" << std::endl;  }
    else                    { glutIdleFunc(idle); std::cout << "middle: off" << std::endl; }
    break;

    case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN) { glutIdleFunc(0);    std::cout << "right: on" << std::endl;  }
    else                    { glutIdleFunc(idle); std::cout << "right: off" << std::endl; }
    break;
  }
}

void resize( int w , int h )
{
  glViewport(0, 0, w, h);

  //glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho( -w / 0.0 , w / 20.0 , -h / 20.0 , h / 20.0 , -3.0 , 3.0 );
  //gluPerspective( 30.0, (double)w / (double)h, 1.0 , 100.0 );
  //gluLookAt( 0.0 , 0.0 , 3.8 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 );
  //monitor.SetWindowSize( w , h );
  //glMatrixMode(GL_MODELVIEW);
}



void OpenGL_init(int *argcp , char **argv)
{
  init();        //初期条件を設定

  glutInit(argcp, argv);

  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

  glutInitWindowSize(monitor.GetWindowSize(Monitor::X),monitor.GetWindowSize(Monitor::Y));
  glutInitWindowPosition( 10 , 100 );//(10,100)
  winid = glutCreateWindow("simulation");
  glutDisplayFunc(display);
  glutReshapeFunc(resize);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glClearColor(1.0, 1.0, 1.0, 1.0);
}

void monitor_init()
{
  double zoom=0.5;
  monitor.SetWindowSize( 800 , 600 );
  //   monitor.SetMode( 0 );
  monitor.SetMovieMode( 1);
  monitor.SetMovieName( "./MovieDir/temp_" );
  monitor.SetZoom(zoom);
  //   monitor.SetGridMode( 0 );
  //   monitor.SetGridWidth( 2.0 );

}

void capture(int *pts)
{
  char filepath[100]; //= "./MovieDir/output.png";
  sprintf(filepath, "./MovieDir/%d.png", *pts);
  png_bytep raw1D;
  png_bytepp raw2D;
  int i;
  int width = glutGet(GLUT_WINDOW_WIDTH);
  int height = glutGet(GLUT_WINDOW_HEIGHT);

  // 構造体確保
  FILE *fp = fopen(filepath, "wb");
  png_structp pp = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  png_infop ip = png_create_info_struct(pp);
  // 書き込み準備
  png_init_io(pp, fp);
  png_set_IHDR(pp, ip, width, height,
               8,                   // 8bit以外にするなら変える
               PNG_COLOR_TYPE_RGBA, // RGBA以外にするなら変える
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // ピクセル領域確保
  raw1D = (png_bytep)malloc(height * png_get_rowbytes(pp, ip));
  raw2D = (png_bytepp)malloc(height * sizeof(png_bytep));
  for (i = 0; i < height; i++)
    raw2D[i] = &raw1D[i * png_get_rowbytes(pp, ip)];
  // 画像のキャプチャ
  glReadBuffer(GL_FRONT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 初期値は4
  glReadPixels(0, 0, width, height,
               GL_RGBA,          // RGBA以外にするなら変える
               GL_UNSIGNED_BYTE, // 8bit以外にするなら変える
               (void *)raw1D);
  // 上下反転
  for (i = 0; i < height / 2; i++)
  {
    png_bytep swp = raw2D[i];
    raw2D[i] = raw2D[height - i - 1];
    raw2D[height - i - 1] = swp;
  }
  // 書き込み
  png_write_info(pp, ip);
  png_write_image(pp, raw2D);
  png_write_end(pp, ip);
  // 開放
  png_destroy_write_struct(&pp, &ip);
  fclose(fp);
  free(raw1D);
  free(raw2D);

  //printf("write out screen capture to '%s'\n", filepath);
}


int main(int argc, char *argv[])
{
  int i, j, i_dim;

  monitor_init();
  std::cout << "monitor init OK" << std::endl;

  OpenGL_init(&argc, argv);

  std::cout << "OpenGL init OK" << std::endl;

  // glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  glutMainLoop();//無限ループ

  return 0;
}
