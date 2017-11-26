#include <Rcpp.h>
#include <math.h>
#include <cmath>

using namespace Rcpp;

const double epsilon = 0.1;

// [[Rcpp::export]]
void motion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, double l, double dt, NumericVector &m, NumericVector &e, double bomb, double energy, double k);

//changing positions and velocities of particles
void euler(NumericVector &q, const NumericVector &v, double dt);

//dealing with wall collisions
void wall_collision(const NumericVector &q, NumericVector &v, double l, NumericVector &m, double k);

//fuse close particles
void fusion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, double l, NumericVector &e, double k);

//explode too big particles
void blust(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, NumericVector &e, double bomb, double energy, double l);

//work out distance
double distance(double x1, double y1, double x2, double y2);

//find first occurence of item in vector
int find(NumericVector &a, int item);

void motion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, double l, double dt, NumericVector &m, NumericVector &e, double bomb, double energy, double k) {
  euler(x,vx,dt);
  euler(y,vy,dt);
  wall_collision(x,vx,l,m,k);
  wall_collision(y,vy,l,m,k);
  fusion(x,y,vx,vy,m,l,e,k);
  blust(x,y,vx,vy,m,e,bomb,energy,l);
  return;
}


void wall_collision(const NumericVector &q, NumericVector &v, double l,
  NumericVector &m, double k) {
  for(int i = 0; i < q.size(); ++i) {
    if(((q[i] < -l + epsilon * m[i] / k) && (v[i] < 0)) ||
                        ((q[i] > l - epsilon * m[i] / k) && (v[i] > 0))) {
      v[i] *= -1;
    }
  }
  return;
}

void fusion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, double l, NumericVector &e, double k) {
  for(int i = 0; i < x.size(); ++i) {
    for(int j = i; j < x.size(); ++j) {
      if(i == j){continue;}
      if(e[i]==0 || e[j]==0){continue;}
      if(distance(x[i],y[i],x[j],y[j]) < epsilon * (m[i] + m[j]) / (k)) {
        int largest = m[i] > m[j] ? i : j;
        int smallest = i + j - largest;
        double mf = m[i] + m[j];
        double vxf = (m[i]*vx[i] + m[j]*vx[j])/mf;
        double vyf = (m[i]*vy[i] + m[j]*vy[j])/mf;
        vx[largest] = vxf;
        vy[largest] = vyf;
        m[largest] = mf;
        m[smallest] = 0;
        vx[smallest] = vy[smallest] = 0;
        x[smallest] = y[smallest] = 2*l;                    //(2l,2l) is point of garbage
        e[smallest] = 0;
      }
    }
  }
  return;
}

void blust(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, NumericVector &e, double bomb, double energy, double l) {
    for(int i = 0; i < x.size(); ++i) {
      if(m[i] >= bomb && e[i]==1) {
        double dalpha = 2 * acos(-1) / m[i];
        double alpha = atan2(vy[i], vx[i]);
        double mi = m[i];
        double dv = energy * m[i];
        double vxi = vx[i];
        double vyi = vy[i];
        double xi = x[i];
        double yi = y[i];
        x[i] = y[i] = 2*l;
        vx[i] = vy[i] = 0;
        m[i] = 0;
        e[i] = 0;
        for(int j = 0; j < mi; ++j) {
          int posj =  find(e, 0);
          vx[posj] = vxi + dv * cos(j * dalpha + alpha);
          vy[posj] = vyi + dv * sin(j * dalpha + alpha);
          x[posj] = xi + mi * epsilon * cos(j * dalpha + alpha) / 2;
          y[posj] = yi + mi * epsilon * sin(j * dalpha + alpha) / 2;
          e[posj] = 1;
          m[posj] = 1;
        }
      }
    }
    return;
  }

double distance(double x1, double y1, double x2, double y2) {
  double dist = pow(x1 - x2, 2) + pow(y1 - y2, 2);
  return pow(dist, 0.5);
}

void euler(NumericVector &q, const NumericVector &v, double dt) {
  for(int i = 0; i < q.size(); ++i) {
    q[i] += v[i] * dt;
  }
  return;
}

int find(NumericVector &a, int item) {
  for(int i = 0; i < a.size(); ++i) {
    if(a[i] == item){return(i);}
  }
  return(-1);
}
