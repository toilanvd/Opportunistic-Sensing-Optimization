#include <bits/stdc++.h>
#include <sys/time.h>

using namespace std;

const double eps = 1e-9;
const double INF = 1e9;

double sqr(double x){return x * x;}
double abso(double x){return (x < 0.0) ? -x : x;}

typedef pair<int, int> ii;
typedef pair<double, double> dd;
typedef pair<double, int> di;
typedef pair<dd, int> intv;

struct busRoute{
    int busId;
    int numTurnPts; //number of turning points
    vector<dd> turnPts; //turning points
};

struct busMap{
    int T1, T2, N;
    vector<busRoute> buses;
    int C;
    vector<ii> criticalSqr; //critical squares
    int unreachableSqr;
};

struct chosenBus{
    int busId;
    vector<ii> coveredSqr; //covered critical squares
};

struct solution{
    int obsSqr; //observable critical squares
    vector<chosenBus> chosenBuses;
};

busMap bmp, originBmp;
double R;
int M, k;

vector<int> criticalSqrCoverCnt;

void readInput(){
    //freopen("Test/10x12_10_0.50.txt", "r", stdin);

    cin >> bmp.T1 >> bmp.T2 >> R >> bmp.N;
    for(int i = 1; i <= bmp.N; i++){
        busRoute bus;
        bus.busId = i;
        cin >> bus.numTurnPts;
        for(int j = 1; j <= bus.numTurnPts; j++){
            double x, y; cin >> x >> y;
            bus.turnPts.push_back(dd(x, y));
        }
        bmp.buses.push_back(bus);
    }
    cin >> bmp.C;
    for(int i = 1; i <= bmp.C; i++){
        int x, y; cin >> x >> y;
        bmp.criticalSqr.push_back(ii(x, y));
        criticalSqrCoverCnt.push_back(0);
    }
    originBmp = bmp;
}

double disPtsToPts(dd p1, dd p2){
    return sqrt(sqr(p1.first - p2.first) + sqr(p1.second - p2.second));
}

double disPtsToSqr(dd p, dd sq){
    //Point inside square => dis = 0
    if(p.first >= sq.first - eps && p.first <= sq.first + 1.0 + eps
       && p.second >= sq.second - eps && p.second <= sq.second + 1.0 + eps){
        return 0.0;
    }
    //Point in x-interval => dis = dis to up/down edge
    if(p.first >= sq.first - eps && p.first <= sq.first + 1.0 + eps){
        return min(abso(p.second - sq.second), abso(p.second - sq.second - 1.0));
    }
    //Point in y-interval => dis = dis to left/right edge
    if(p.second >= sq.second - eps && p.second <= sq.second + 1.0 + eps){
        return min(abso(p.first - sq.first), abso(p.first - sq.first - 1.0));
    }
    //Else, dis = min dis to 4 vertices
    double dis = disPtsToPts(p, dd(sq.first, sq.second));
    dis = min(dis, disPtsToPts(p, dd(sq.first + 1.0, sq.second)));
    dis = min(dis, disPtsToPts(p, dd(sq.first, sq.second + 1.0)));
    dis = min(dis, disPtsToPts(p, dd(sq.first + 1.0, sq.second + 1.0)));
    return dis;
}

bool ptsInRec(dd p, dd q1, dd q2){
    if(((q1.first - eps <= p.first && p.first <= q2.first + eps) || (q1.first + eps >= p.first && p.first >= q2.first - eps))
       && ((q1.second - eps <= p.second && p.second <= q2.second + eps) || (q1.second + eps >= p.second && p.second >= q2.second - eps))){
        return true;
    }
    return false;
}

bool ptsInQuarter(dd p, dd O, int dir){
    if(dir == 0 && p.first <= O.first + eps && p.second + eps >= O.second) return true;
    if(dir == 1 && p.first + eps >= O.first && p.second + eps >= O.second) return true;
    if(dir == 2 && p.first + eps >= O.first && p.second <= O.second + eps) return true;
    if(dir == 3 && p.first <= O.first + eps && p.second <= O.second + eps) return true;
    return false;
}

vector<double> intsSegAndSeg(dd p1, dd p2, dd q1, dd q2){
    vector<double> res;
    //ax + by = c
    double a = p2.second - p1.second;
    double b = -(p2.first - p1.first);
    double c = (p2.second - p1.second) * p1.first - (p2.first - p1.first) * p1.second;
    if(q1.first == q2.first){ //vertical segment: x = d
        double d = q1.first;
        if(abso(b) <= eps){ //parallel lines
            if(abso(c/a - d) <= eps){ //same lines
                double l = max(min(p1.second, p2.second), q1.second);
                double r = min(max(p1.second, p2.second), q2.second);
                if(l <= r + eps){
                    res.push_back(disPtsToPts(p1, dd(d, l)));
                    res.push_back(disPtsToPts(p1, dd(d, r)));
                }
            }
        }
        else{ //non-parallel lines
            double y = (c - a * d) / b;
            if(y >= q1.second - eps && y <= q2.second + eps
               && ptsInRec(dd(d, y), p1, p2)){
                res.push_back(disPtsToPts(p1, dd(d, y)));
            }
        }
    }
    else{ //horizontal segment: y = d
        double d = q1.second;
        if(abso(a) <= eps){ //parallel lines
            if(abso(c/b - d) <= eps){ //same lines
                double l = max(min(p1.first, p2.first), q1.first);
                double r = min(max(p1.first, p2.first), q2.first);
                if(l <= r + eps){
                    res.push_back(disPtsToPts(p1, dd(l, d)));
                    res.push_back(disPtsToPts(p1, dd(r, d)));
                }
            }
        }
        else{ //non-parallel lines
            double x = (c - b * d) / a;
            if(x >= q1.first - eps && x <= q2.first + eps
               && ptsInRec(dd(x, d), p1, p2)){
                res.push_back(disPtsToPts(p1, dd(x, d)));
            }
        }
    }

    return res;
}

vector<double> intsSegAndQuarterArc(dd p1, dd p2, dd O, int dir){
    //dir: the quarter which the arc belongs to (0: north-west, 1: north-east, 2: south-east, 3: south-west)
    vector<double> res;
    //ax + by + c = 0
    double a = p2.second - p1.second;
    double b = -(p2.first - p1.first);
    double c = (p2.first - p1.first) * (p1.second - O.second) - (p2.second - p1.second) * (p1.first - O.first);

    double x0 = -a*c/(a*a+b*b), y0 = -b*c/(a*a+b*b);
    if (abs (c*c - R*R*(a*a+b*b)) <= eps) { //tangent
        x0 += O.first;
        y0 += O.second;
        if(ptsInRec(dd(x0, y0), p1, p2) && ptsInQuarter(dd(x0, y0), O, dir)){
            res.push_back(disPtsToPts(p1, dd(x0, y0)));
        }
        //cout << x0 << ' ' << y0 << '\n';
    }
    else if(c*c + eps <= R*R*(a*a+b*b)){ //2 points of intersection
        double d = R*R - c*c/(a*a+b*b);
        double mult = sqrt (d / (a*a+b*b));
        double ax, ay, bx, by;
        ax = x0 + b * mult + O.first;
        bx = x0 - b * mult + O.first;
        ay = y0 - a * mult + O.second;
        by = y0 + a * mult + O.second;
        if(ptsInRec(dd(ax, ay), p1, p2) && ptsInQuarter(dd(ax, ay), O, dir)){
            res.push_back(disPtsToPts(p1, dd(ax, ay)));
        }
        if(ptsInRec(dd(bx, by), p1, p2) && ptsInQuarter(dd(bx, by), O, dir)){
            res.push_back(disPtsToPts(p1, dd(bx, by)));
        }
        //cout << ax << ' ' << ay << '\n' << bx << ' ' << by << '\n';
    }

    return res;
}

//Calculate intersections between segment and extended square
vector<double> intsSegAndExtSqr(dd p1, dd p2, dd sq, double lenFromStart){
    vector<double> res;
    vector<double> temp;
    temp = intsSegAndSeg(p1, p2, dd(sq.first, sq.second - R), dd(sq.first + 1.0, sq.second - R));
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    temp = intsSegAndSeg(p1, p2, dd(sq.first - R, sq.second), dd(sq.first - R, sq.second + 1.0));
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    temp = intsSegAndSeg(p1, p2, dd(sq.first, sq.second + 1.0 + R), dd(sq.first + 1.0, sq.second + 1.0 + R));
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    temp = intsSegAndSeg(p1, p2, dd(sq.first + 1.0 + R, sq.second), dd(sq.first + 1.0 + R, sq.second + 1.0));
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    temp = intsSegAndQuarterArc(p1, p2, dd(sq.first, sq.second + 1.0), 0);
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    temp = intsSegAndQuarterArc(p1, p2, dd(sq.first + 1.0, sq.second + 1.0), 1);
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    temp = intsSegAndQuarterArc(p1, p2, dd(sq.first + 1.0, sq.second), 2);
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    temp = intsSegAndQuarterArc(p1, p2, dd(sq.first, sq.second), 3);
    while(!temp.empty()){res.push_back(lenFromStart + temp.back()); temp.pop_back();}
    return res;
}

int g[5005][5005];
int f[5005][2505];

chosenBus optimalSet(int busIdInCurMap){
    chosenBus cb;
    busRoute Y = bmp.buses[busIdInCurMap];
    cb.busId = Y.busId;

    //Calculate the observable intervals
    vector<intv> obsIntv;
    vector<double> criticalPts;
    criticalPts.push_back(-INF);
    for(int i = 0; i < bmp.C; i++){
        dd X = dd((double)bmp.criticalSqr[i].first, (double)bmp.criticalSqr[i].second);
        vector<double> L;
        double lenFromStart = 0.0;
        //If starting point of Y can observe X
        if(disPtsToSqr(Y.turnPts[0], X) <= R + eps){
            L.push_back(0.0);
        }
        //Loop through all intervals in Y
        for(int j = 0; j < Y.numTurnPts - 1; j++){
            vector<double> tempIntersections = intsSegAndExtSqr(Y.turnPts[j], Y.turnPts[j+1], X, lenFromStart);
            while(!tempIntersections.empty()){
                L.push_back(tempIntersections.back());
                tempIntersections.pop_back();
            }
            lenFromStart += disPtsToPts(Y.turnPts[j], Y.turnPts[j+1]);
        }
        //If ending point of Y can observe X
        if(Y.numTurnPts > 1 && disPtsToSqr(Y.turnPts[Y.numTurnPts - 1], X) <= R + eps){
            L.push_back(lenFromStart);
        }
        //Get observable interval
        if(!L.empty()){
            criticalSqrCoverCnt[i]++;
            //if(i+1==10) cout << i+1 << " - " << busIdInCurMap+1 << endl;
            sort(L.begin(), L.end());
            vector<double> tempL;
            tempL.push_back(L[0]);
            for(int j = 1; j < (int)L.size(); j++){
                if(abso(L[j] - L[j - 1]) > eps){
                    tempL.push_back(L[j]);
                }
            }
            L = tempL;
            if((int)L.size() == 1) L.push_back(L[0]);
            obsIntv.push_back(intv(dd(L[0], L[(int)L.size()-1]), i));
            criticalPts.push_back(L[0]);
            criticalPts.push_back(L[(int)L.size()-1]);
        }
    }

    //Sort 2D critical points in increasing order
    sort(obsIntv.begin(), obsIntv.end());
    sort(criticalPts.begin(), criticalPts.end());
    int D = (int)obsIntv.size();

    //Calculate function g(i, j) for all 0 <= i < j <= 2D
    for(int j = 0; j <= 2 * D; j++){
        int cnt = 0;
        for(int i = 0; i < D; i++){
            if(criticalPts[j] + eps < obsIntv[i].first.first) break;
            if(criticalPts[j] >= obsIntv[i].first.first - eps
               && criticalPts[j] <= obsIntv[i].first.second + eps){
                cnt++;
            }
        }

        int t = 0;
        for(int i = 0; i < j; i++){
            while(t < D && criticalPts[i] >= obsIntv[t].first.first - eps){
                if(criticalPts[j] <= obsIntv[t].first.second + eps){
                    cnt--;
                }
                t++;
            }
            g[i][j] = cnt;
        }
    }

    //Calculate function f(i, j) for all 0 <= i < j <= 2D
    for(int j = 0; j <= 2 * k; j++){
        for(int i = 0; i <= 2 * D; i++){
            f[i][j] = 0;
            if(j > 0){
                f[i][j] = f[i][j - 1];
                for(int t = i + 1; t <= 2 * D; t++){
                    int v = f[t][j - 1] + g[i][t];
                    if(f[i][j] < v){
                        f[i][j] = v;
                    }
                }
            }
        }
    }

    //Trace result from function f
    vector<double> turnedOnPos;
    int x = 0, y = 2 * k;
    while(y > 0){
        if(f[x][y] == f[x][y - 1]) y--;
        else{
            for(int z = x + 1; z <= 2 * D; z++){
                if(f[x][y] == f[z][y - 1] + g[x][z]){
                    turnedOnPos.push_back(criticalPts[z]);
                    x = z; y--;
                    break;
                }
            }
        }
    }
    for(int i = 0; i < D; i++){
        for(int j = 0; j < (int)turnedOnPos.size(); j++){
            if(obsIntv[i].first.first - eps <= turnedOnPos[j]
               && obsIntv[i].first.second + eps >= turnedOnPos[j]){
                cb.coveredSqr.push_back(bmp.criticalSqr[obsIntv[i].second]);
                break;
            }
        }
    }

    return cb;
}

vector<ii> fullSet(int busIdInCurMap){
    vector<ii> res;
    busRoute Y = bmp.buses[busIdInCurMap];

    //Calculate the observable intervals
    for(int i = 0; i < bmp.C; i++){
        dd X = dd((double)bmp.criticalSqr[i].first, (double)bmp.criticalSqr[i].second);
        vector<double> L;
        double lenFromStart = 0.0;
        //If starting point of Y can observe X
        if(disPtsToSqr(Y.turnPts[0], X) <= R + eps){
            L.push_back(0.0);
        }
        //Loop through all intervals in Y
        for(int j = 0; j < Y.numTurnPts - 1; j++){
            vector<double> tempIntersections = intsSegAndExtSqr(Y.turnPts[j], Y.turnPts[j+1], X, lenFromStart);
            while(!tempIntersections.empty()){
                L.push_back(tempIntersections.back());
                tempIntersections.pop_back();
            }
            lenFromStart += disPtsToPts(Y.turnPts[j], Y.turnPts[j+1]);
        }
        //If ending point of Y can observe X
        if(Y.numTurnPts > 1 && disPtsToSqr(Y.turnPts[Y.numTurnPts - 1], X) <= R + eps){
            L.push_back(lenFromStart);
        }
        //Get observable interval
        if(!L.empty()) res.push_back(bmp.criticalSqr[i]);
    }

    return res;
}

void calculateUnreachableSqr(){
    fill(criticalSqrCoverCnt.begin(), criticalSqrCoverCnt.end(), 0);

    for(int i = 0; i < bmp.N; i++) optimalSet(i);

    originBmp.unreachableSqr = 0;
    for(int i = 0; i < bmp.C; i++) originBmp.unreachableSqr += (criticalSqrCoverCnt[i] == 0) ? 1 : 0;
}

int upperboundOPT(){
    vector<int> cvSqr;
    for(int i = 0; i < bmp.N; i++){
        chosenBus X = optimalSet(i);
        cvSqr.push_back((int)X.coveredSqr.size());
    }
    sort(cvSqr.begin(), cvSqr.end());
    int upbOPT = 0;
    for(int i = bmp.N - 1; i >= bmp.N - M; i--){
        upbOPT += cvSqr[i];
    }
    return min(upbOPT, bmp.C - bmp.unreachableSqr);
}

int maxCoveredSqr(){
    fill(criticalSqrCoverCnt.begin(), criticalSqrCoverCnt.end(), 0);

    for(int i = 0; i < bmp.N; i++) optimalSet(i);

    int res = 0;
    for(int i = 0; i < bmp.C; i++) res = max(res, criticalSqrCoverCnt[i]);

    return res;
}

double avgCoveredSqr(){
    int res = 0;
    for(int i = 0; i < bmp.C; i++) res += criticalSqrCoverCnt[i];

    return (double)res / (double)bmp.C;
}

double avgCommonSqr(){
    vector<vector<ii> > tempVec;
    for(int i = 0; i < bmp.N; i++) tempVec.push_back(fullSet(i));
    int res = 0;
    for(int i = 0; i < bmp.N; i++){
        for(int j = 0; j < i; j++){
            for(int t1 = 0; t1 < (int)tempVec[i].size(); t1++){
                for(int t2 = 0; t2 < (int)tempVec[j].size(); t2++){
                    if(tempVec[i][t1] == tempVec[j][t2]) res++;
                }
            }
        }
    }

    return (double)res / (double)((bmp.N*(bmp.N-1))/2);
}

solution solveSOBP(){
    solution S;
    S.obsSqr = 0;
    for(int i = 1; i <= M; i++){
        chosenBus X;
        int chosenIdInCurMap;
        for(int j = 0; j < bmp.N; j++){
            chosenBus Y = optimalSet(j);
            if((int)Y.coveredSqr.size() >= (int)X.coveredSqr.size()){
                X = Y;
                chosenIdInCurMap = j;
            }
        }
        S.chosenBuses.push_back(X);
        S.obsSqr += (int)X.coveredSqr.size();
        bmp.buses.erase(bmp.buses.begin() + chosenIdInCurMap);
        bmp.N--;
        for(int j = 0; j < (int)X.coveredSqr.size(); j++){
            for(int t = 0; t < (int)bmp.criticalSqr.size(); t++){
                if(X.coveredSqr[j] == bmp.criticalSqr[t]){
                    bmp.criticalSqr.erase(bmp.criticalSqr.begin() + t);
                    break;
                }
            }
        }
        bmp.C -= (int)X.coveredSqr.size();
    }
    return S;
}

void printResult(solution res){
    cout << "Covered critical squares: " << res.obsSqr << endl;
    cout << "Chosen buses:" << endl;
    for(int i = 0; i < M; i++){
        cout << res.chosenBuses[i].busId << ":";
        for(int j = 0; j < (int)res.chosenBuses[i].coveredSqr.size(); j++){
            cout << " (" << res.chosenBuses[i].coveredSqr[j].first;
            cout << ", " << res.chosenBuses[i].coveredSqr[j].second << "),";
        }
        cout << endl;
    }
}

// Check all cases of M and k
void fullCheck(){
    cout.precision(3);
    ///* Time calculation
    double mapRunningTime = 0.0;
    //*/
    /* Early termination
    int tempCnt = 0;
    //*/
    int test_cnt = 0;
    double minEfficiency = 100.0;
    for(k = 1; k <= originBmp.C; k++){
        for(M = 1; M <= originBmp.N; M++){
            test_cnt++;
            ///* Time calculation
            struct timeval startTestTime, endTestTime;
            gettimeofday(&startTestTime, NULL);
            //*/
            bmp = originBmp;
            int upbOPT = upperboundOPT();
            solution res = solveSOBP();
            double eff = 100.0 * (double)res.obsSqr / (double)upbOPT;
            /* Early termination
            if(eff + eps > minEfficiency) tempCnt++; else tempCnt = 0;
            if(tempCnt >= 280) break;
            //*/
            minEfficiency = min(minEfficiency, eff);
            ///* Time calculation
            gettimeofday(&endTestTime, NULL);
            double testRunningTime = double(endTestTime.tv_sec - startTestTime.tv_sec) + double(endTestTime.tv_usec - startTestTime.tv_usec) * 1e-6;
            mapRunningTime += testRunningTime;
            //*/
            ///* Print test result
            //printResult(res);
            cout << "T1 = " << originBmp.T1 << ", T2 = " << originBmp.T2 << ", C = " << originBmp.C << ", R = " << fixed << R << ", ";
            cout << "M = " << M << ", k = " << k << ": upbOPT = " << upbOPT << ", res = " << res.obsSqr;
            cout << " | Eff: " << fixed << eff << ", minEff: " << fixed << minEfficiency;
            if(abso(eff - minEfficiency) <= eps) cout << " (***)";
            cout << ", time: " << fixed << testRunningTime << " seconds";
            cout << endl;
            //*/
            /* Early termination
            if(upbOPT == originBmp.C - originBmp.unreachableSqr) break;
            //*/
        }
        /* Early termination
        if(tempCnt >= 280) break;
        //*/
    }

    cout << "Minimum efficiency: " << fixed << minEfficiency << "%" << endl;
    ///* Time calculation
    cout << "Time taken: " << fixed << mapRunningTime << " seconds" << endl;
    cout << "Avg time per (M, k): " << fixed << mapRunningTime / (double)test_cnt << " seconds" << endl;
    //*/
}

int main(){
    readInput();
    calculateUnreachableSqr();
    /*
    cout << "max covered square = " << maxCoveredSqr();
    cout << ", average covered = " << avgCoveredSqr();
    cout << ", average common = " << (double)avgCommonSqr() / ((double)bmp.C / (double)(bmp.T1 * bmp.T2)) << endl;
    */
    fullCheck();

    /*Test
    M = 8; k = 1;
    bmp = originBmp;
    int upbOPT = upperboundOPT();
    cout << "Upper bound of OPT = " << upbOPT << endl;
    solution res = solveSOBP();
    printResult(res);
    double eff = 100.0 * (double)res.obsSqr / (double)upbOPT;
    cout << "Efficiency >= " << eff << endl;
    //*/

    return 0;
}
