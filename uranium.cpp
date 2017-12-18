#include <Rcpp.h>
#include <math.h>
#include <cmath>

using namespace Rcpp;

const double epsilon = 0.07;

// [[Rcpp::export]]
void motion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, double l, double dt, NumericVector &m, LogicalVector &e, double bomb, double energy, double k);

//changing positions and velocities of particles
void euler(NumericVector &q, const NumericVector &v, double dt);

//dealing with wall collisions
void wall_collision(const NumericVector &q, NumericVector &v, double l, NumericVector &m, double k);

//fuse close particles
void fusion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, double l, LogicalVector &e, double k);

//explode too big particles
void fission(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, LogicalVector &e, double bomb, double energy, double l);

//work out distance
double distance(double x1, double y1, double x2, double y2);

//find first occurence of item in vector
int find_false(const LogicalVector &a);

void motion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, double l, double dt, NumericVector &m, LogicalVector &e, double bomb, double energy, double k) {
  euler(x,vx,dt);
  euler(y,vy,dt);
  wall_collision(x,vx,l,m,k);
  wall_collision(y,vy,l,m,k);
  fusion(x,y,vx,vy,m,l,e,k);
  fission(x,y,vx,vy,m,e,bomb,energy,l);
}


void wall_collision(const NumericVector &q, NumericVector &v, double l,
  NumericVector &m, double k) {
  for(int i = 0; i < q.size(); ++i) {
    if(((q[i] < -l + epsilon * m[i] / k) && (v[i] < 0)) ||
                        ((q[i] > l - epsilon * m[i] / k) && (v[i] > 0))) {
      v[i] *= -1;
    }
  }
}

void fusion(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, double l, LogicalVector &e, double k) {
  for(int i = 0; i < x.size(); ++i) {
    for(int j = i; j < x.size(); ++j) {
      if(i == j){continue;}
      if(!e[i] || !e[j]){continue;}
      if(distance(x[i],y[i],x[j],y[j]) < epsilon * (m[i] + m[j]) / k) {
        int lt = m[i] > m[j] ? i : j;
        int sl = i + j - lt;
        double mf = m[i] + m[j];
        double vxf = (m[i]*vx[i] + m[j]*vx[j])/mf, vyf = (m[i]*vy[i] + m[j]*vy[j])/mf;
        vx[lt] = vxf;
        vy[lt] = vyf;
        m[lt] = mf;
        e[sl] = false;
      }
    }
  }
}

void fission(NumericVector &x, NumericVector &y, NumericVector &vx,
  NumericVector &vy, NumericVector &m, LogicalVector &e, double bomb, double energy, double l) {
    for(int i = 0; i < x.size(); ++i) {
      if(m[i] >= bomb && e[i]) {
        double delta = 2 * acos(-1) / m[i], alpha = atan2(vy[i], vx[i]);
        double vxi = vx[i], vyi = vy[i], xi = x[i], yi = y[i], dv = energy * m[i], mi = m[i];
        e[i] = false;
        for(int j = 0; j < mi; ++j) {
          int posj =  find_false(e);
          vx[posj] = vxi + dv * cos(j * delta + alpha);
          vy[posj] = vyi + dv * sin(j * delta + alpha);
          x[posj] = xi + mi * epsilon * cos(j * delta + alpha) / 2;
          y[posj] = yi + mi * epsilon * sin(j * delta + alpha) / 2;
          e[posj] = true;
          m[posj] = 1;
        }
      }
    }
  }

double distance(double x1, double y1, double x2, double y2) {
  double dist = pow(x1 - x2, 2) + pow(y1 - y2, 2);
  return pow(dist, 0.5);
}

void euler(NumericVector &q, const NumericVector &v, double dt) {
  for(int i = 0; i < q.size(); ++i) {
    q[i] += v[i] * dt;
  }
}

int find_false(const LogicalVector &a) {
  for(int i = 0; i < a.size(); ++i) {
    if(!a[i]){return(i);}
  }
  return(-1);
}
