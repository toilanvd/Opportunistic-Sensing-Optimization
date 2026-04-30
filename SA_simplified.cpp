#include <bits/stdc++.h>
#include <sys/time.h>

using namespace std;

const double EPS = 1e-9;
const double INF = 1e9;

const double T_MAX = 100.0;
const double T_DEC = 0.1;
const double T_MIN = 0.1;
const int LOOPTIME = 200;
const bool USE_GREEDY_INITIAL_SOLUTION = false;

double sqr(double x){return x * x;}
double abso(double x){return (x < 0.0) ? -x : x;}

typedef pair<int, int> ii;
typedef pair<double, double> dd;
typedef pair<double, int> di;
typedef pair<dd, int> intv;
typedef pair<ii, int> iii;

struct busRoute{
    int busId;
    int numTurnPts; //number of turning points
    vector<dd> turnPts; //turning points
    vector<vector<ii> > obsCritSqr; // * observable squares by a critical point
};

struct busMap{
    int T1, T2, N;
    vector<busRoute> buses;
    int C;
    vector<ii> criticalSqr; //critical squares
    int unreachableSqrCnt;
};

struct chosenBus{
    int busId;
    vector<ii> coveredSqr; //covered critical squares
    vector<int> critPtsSet; // * set of index of chosen critical points
};

struct solution{
    int obsSqrCnt; //observable critical squares
    vector<chosenBus> chosenBuses;
};

busMap bmp, originBmp;
double R;
int M, k;

vector<int> criticalSqrCoverCnt;

void readInput(){
    //freopen("Test_simplified/10x12_10_0.50.txt", "r", stdin);

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
    if(p.first >= sq.first - EPS && p.first <= sq.first + 1.0 + EPS
       && p.second >= sq.second - EPS && p.second <= sq.second + 1.0 + EPS){
        return 0.0;
    }
    //Point in x-interval => dis = dis to up/down edge
    if(p.first >= sq.first - EPS && p.first <= sq.first + 1.0 + EPS){
        return min(abso(p.second - sq.second), abso(p.second - sq.second - 1.0));
    }
    //Point in y-interval => dis = dis to left/right edge
    if(p.second >= sq.second - EPS && p.second <= sq.second + 1.0 + EPS){
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
    if(((q1.first - EPS <= p.first && p.first <= q2.first + EPS) || (q1.first + EPS >= p.first && p.first >= q2.first - EPS))
       && ((q1.second - EPS <= p.second && p.second <= q2.second + EPS) || (q1.second + EPS >= p.second && p.second >= q2.second - EPS))){
        return true;
    }
    return false;
}

bool ptsInQuarter(dd p, dd O, int dir){
    if(dir == 0 && p.first <= O.first + EPS && p.second + EPS >= O.second) return true;
    if(dir == 1 && p.first + EPS >= O.first && p.second + EPS >= O.second) return true;
    if(dir == 2 && p.first + EPS >= O.first && p.second <= O.second + EPS) return true;
    if(dir == 3 && p.first <= O.first + EPS && p.second <= O.second + EPS) return true;
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
        if(abso(b) <= EPS){ //parallel lines
            if(abso(c/a - d) <= EPS){ //same lines
                double l = max(min(p1.second, p2.second), q1.second);
                double r = min(max(p1.second, p2.second), q2.second);
                if(l <= r + EPS){
                    res.push_back(disPtsToPts(p1, dd(d, l)));
                    res.push_back(disPtsToPts(p1, dd(d, r)));
                }
            }
        }
        else{ //non-parallel lines
            double y = (c - a * d) / b;
            if(y >= q1.second - EPS && y <= q2.second + EPS
               && ptsInRec(dd(d, y), p1, p2)){
                res.push_back(disPtsToPts(p1, dd(d, y)));
            }
        }
    }
    else{ //horizontal segment: y = d
        double d = q1.second;
        if(abso(a) <= EPS){ //parallel lines
            if(abso(c/b - d) <= EPS){ //same lines
                double l = max(min(p1.first, p2.first), q1.first);
                double r = min(max(p1.first, p2.first), q2.first);
                if(l <= r + EPS){
                    res.push_back(disPtsToPts(p1, dd(l, d)));
                    res.push_back(disPtsToPts(p1, dd(r, d)));
                }
            }
        }
        else{ //non-parallel lines
            double x = (c - b * d) / a;
            if(x >= q1.first - EPS && x <= q2.first + EPS
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
    if (abs (c*c - R*R*(a*a+b*b)) <= EPS) { //tangent
        x0 += O.first;
        y0 += O.second;
        if(ptsInRec(dd(x0, y0), p1, p2) && ptsInQuarter(dd(x0, y0), O, dir)){
            res.push_back(disPtsToPts(p1, dd(x0, y0)));
        }
        //cout << x0 << ' ' << y0 << '\n';
    }
    else if(c*c + EPS <= R*R*(a*a+b*b)){ //2 points of intersection
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
        if(disPtsToSqr(Y.turnPts[0], X) <= R + EPS){
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
        if(Y.numTurnPts > 1 && disPtsToSqr(Y.turnPts[Y.numTurnPts - 1], X) <= R + EPS){
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
                if(abso(L[j] - L[j - 1]) > EPS){
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
            if(criticalPts[j] + EPS < obsIntv[i].first.first) break;
            if(criticalPts[j] >= obsIntv[i].first.first - EPS
               && criticalPts[j] <= obsIntv[i].first.second + EPS){
                cnt++;
            }
        }

        int t = 0;
        for(int i = 0; i < j; i++){
            while(t < D && criticalPts[i] >= obsIntv[t].first.first - EPS){
                if(criticalPts[j] <= obsIntv[t].first.second + EPS){
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
            if(obsIntv[i].first.first - EPS <= turnedOnPos[j]
               && obsIntv[i].first.second + EPS >= turnedOnPos[j]){
                cb.coveredSqr.push_back(bmp.criticalSqr[obsIntv[i].second]);
                break;
            }
        }
    }

    return cb;
}

void buildObsCritSqrForBus(int busIdInCurMap){
    busRoute Y = originBmp.buses[busIdInCurMap];

    //Calculate the observable intervals (critical points are endpoints of these intervals)
    vector<intv> obsIntv;
    for(int i = 0; i < originBmp.C; i++){
        dd X = dd((double)originBmp.criticalSqr[i].first, (double)originBmp.criticalSqr[i].second);
        vector<double> L;
        double lenFromStart = 0.0;
        //If starting point of Y can observe X
        if(disPtsToSqr(Y.turnPts[0], X) <= R + EPS){
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
        if(Y.numTurnPts > 1 && disPtsToSqr(Y.turnPts[Y.numTurnPts - 1], X) <= R + EPS){
            L.push_back(lenFromStart);
        }
        //Get observable interval
        if(!L.empty()){
            sort(L.begin(), L.end());
            vector<double> tempL;
            tempL.push_back(L[0]);
            for(int j = 1; j < (int)L.size(); j++){
                if(abso(L[j] - L[j - 1]) > EPS){
                    tempL.push_back(L[j]);
                }
            }
            L = tempL;
            if((int)L.size() == 1) L.push_back(L[0]);
            obsIntv.push_back(intv(dd(L[0], L[(int)L.size()-1]), i));
        }
    }

    // Calculate observable critical square set from each critical point
    for(int i = 0; i < (int)obsIntv.size(); i++){
        vector<ii> tempObsCritSqr[2];
        for(int j = 0; j < (int)obsIntv.size(); j++){
            if(obsIntv[i].first.first + EPS >= obsIntv[j].first.first && obsIntv[i].first.first <= obsIntv[j].first.second + EPS){
                tempObsCritSqr[0].push_back(originBmp.criticalSqr[obsIntv[j].second]);
            }
            if(obsIntv[i].first.second + EPS >= obsIntv[j].first.first && obsIntv[i].first.second <= obsIntv[j].first.second + EPS){
                tempObsCritSqr[1].push_back(originBmp.criticalSqr[obsIntv[j].second]);
            }
        }
        originBmp.buses[busIdInCurMap].obsCritSqr.push_back(tempObsCritSqr[0]);
        originBmp.buses[busIdInCurMap].obsCritSqr.push_back(tempObsCritSqr[1]);
    }
}

void buildObsCritSqr(){
    for(int i = 0; i < originBmp.N; i++){
        buildObsCritSqrForBus(i);
    }
}

void calCoveredSqr(chosenBus* X){
    set<ii> obsSqr;
    for(int i = 0; i < originBmp.N; i++){
        if(originBmp.buses[i].busId == (*X).busId){
            for(int j = 0; j < (int)(*X).critPtsSet.size(); j++){
                for(int t = 0; t < (int)originBmp.buses[i].obsCritSqr[(*X).critPtsSet[j]].size(); t++){
                    obsSqr.insert(originBmp.buses[i].obsCritSqr[(*X).critPtsSet[j]][t]);
                }
            }
            break;
        }
    }
    vector<ii> tempCoveredSqr;
    for(set<ii> :: iterator it = obsSqr.begin(); it != obsSqr.end(); ++it){
        tempCoveredSqr.push_back((*it));
    }
    (*X).coveredSqr = tempCoveredSqr;
}

void calSolutionStats(solution* sol){
    // Must calculated coveredSqr for each bus route beforehand
    set<ii> obsSqr;
    for(int i = 0; i < (int)(*sol).chosenBuses.size(); i++){
        for(int j = 0; j < (int)(*sol).chosenBuses[i].coveredSqr.size(); j++){
            obsSqr.insert((*sol).chosenBuses[i].coveredSqr[j]);
        }
    }
    (*sol).obsSqrCnt = (int)obsSqr.size();
}

void calculateUnreachableSqr(){
    fill(criticalSqrCoverCnt.begin(), criticalSqrCoverCnt.end(), 0);

    for(int i = 0; i < originBmp.N; i++) optimalSet(i);

    originBmp.unreachableSqrCnt = 0;
    for(int i = 0; i < originBmp.C; i++) originBmp.unreachableSqrCnt += (criticalSqrCoverCnt[i] == 0) ? 1 : 0;
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
    return min(upbOPT, bmp.C - bmp.unreachableSqrCnt);
}

chosenBus createRandomizedChosenBus(int busId){
    chosenBus res;
    res.busId = busId;
    int index = -1;
    for(int i = 0; i < originBmp.N; i++){
        if(originBmp.buses[i].busId == busId){
            index = i;
            break;
        }
    }
    vector<bool> pickedCritPts((int)originBmp.buses[index].obsCritSqr.size(), false);
    for(int i = 1; i <= min(2 * k, (int)originBmp.buses[index].obsCritSqr.size()); i++){
        int cnt = rand() % ((int)originBmp.buses[index].obsCritSqr.size() - i + 1) + 1;
        for(int j = 0; j < (int)pickedCritPts.size(); j++){
            if(!pickedCritPts[j]){
                cnt--;
                if(cnt == 0){
                    pickedCritPts[j] = true;
                    break;
                }
            }
        }
    }
    for(int i = 0; i < (int)pickedCritPts.size(); i++) if(pickedCritPts[i]) res.critPtsSet.push_back(i);
    calCoveredSqr(&res);
    return res;
}

solution createRandomizedSolution(){
    solution S;
    S.obsSqrCnt = 0;
    vector<bool> pickedBuses(originBmp.N, false);
    for(int t = 0; t < M; t++){
        int cnt = rand() % (originBmp.N - t) + 1;
        for(int j = 0; j < (int)pickedBuses.size(); j++){
            if(!pickedBuses[j]){
                cnt--;
                if(cnt == 0){
                    pickedBuses[j] = true;
                    break;
                }
            }
        }
    }
    for(int j = 0; j < (int)pickedBuses.size(); j++) if(pickedBuses[j]) S.chosenBuses.push_back(createRandomizedChosenBus(originBmp.buses[j].busId));
    calSolutionStats(&S);
    return S;
}

void printSolution(solution res){
    cout << "Covered critical squares: " << res.obsSqrCnt << endl;
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

bool cmpSolution(solution A, solution B){
    if(A.obsSqrCnt > B.obsSqrCnt) return true;
    else return false;
}

busMap makeNewBmp(solution S){
    busMap res = originBmp;

    for(int i = 0; i < (int)S.chosenBuses.size(); i++){
        res.N--;
        for(int j = 0; j < (int)res.buses.size(); j++){
            if(res.buses[j].busId == S.chosenBuses[i].busId){
                res.buses.erase(res.buses.begin() + j);
                break;
            }
        }
        for(int j = 0; j < (int)S.chosenBuses[i].coveredSqr.size(); j++){
            for(int t = 0; t < (int)res.criticalSqr.size(); t++){
                if(S.chosenBuses[i].coveredSqr[j] == res.criticalSqr[t]){
                    res.criticalSqr.erase(res.criticalSqr.begin() + t);
                    break;
                }
            }
        }
        res.C = (int)res.criticalSqr.size();
    }

    return res;
}

chosenBus createGreedyChosenBus(int busId){
    for(int i = 0; i < (int)bmp.buses.size(); i++){
        if(bmp.buses[i].busId == busId){
            return optimalSet(i);
        }
    }
}

int currentCoverSqrCnt(chosenBus X){
    int res = 0;
    for(int i = 0; i < (int)X.coveredSqr.size(); i++){
        for(int j = 0; j < (int)bmp.criticalSqr.size(); j++){
            if(X.coveredSqr[i] == bmp.criticalSqr[j]){
                res++;
                break;
            }
        }
    }
    return res;
}

solution createGreedyInitialSolution(){
    solution S;
    S.obsSqrCnt = 0;
    for(int i = 1; i <= M; i++){
        chosenBus X;
        int chosenIdInCurMap;
        for(int j = 0; j < bmp.N; j++){
            chosenBus Y = optimalSet(j);
            if(currentCoverSqrCnt(Y) >= currentCoverSqrCnt(X)){
                X = Y;
                chosenIdInCurMap = j;
            }
        }
        S.chosenBuses.push_back(X);
        S.obsSqrCnt += (int)X.coveredSqr.size();
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
        bmp.C = (int)bmp.criticalSqr.size();
    }
    calSolutionStats(&S);
    bmp = originBmp;
    return S;
}

solution solveSOBP(){
    /// SA greedyK simplified algorithm

    solution S = USE_GREEDY_INITIAL_SOLUTION ? createGreedyInitialSolution() : createRandomizedSolution();
    solution res = S;
    double T = T_MAX;
    while(T + EPS >= T_MIN){
        for(int timer = 1; timer <= LOOPTIME; timer++){
            solution newS = S;
            vector<bool> usedBus(originBmp.N + 1, false);
            for(int j = 0; j < M; j++) usedBus[S.chosenBuses[j].busId] = true;
            int busToRemove = rand() % M;
            usedBus[S.chosenBuses[busToRemove].busId] = false;
            newS.chosenBuses.erase(newS.chosenBuses.begin() + busToRemove);
            int cntBusToAdd = rand() % (originBmp.N - M + 1) + 1, busToAdd = -1;
            for(int j = 0; j < originBmp.N; j++){
                if(!usedBus[originBmp.buses[j].busId]){
                    cntBusToAdd--;
                    if(cntBusToAdd == 0){
                        busToAdd = originBmp.buses[j].busId;
                        break;
                    }
                }
            }
            bmp = makeNewBmp(newS);
            newS.chosenBuses.push_back(createGreedyChosenBus(busToAdd));
            calSolutionStats(&newS);
            if(cmpSolution(newS, res)) res = newS;
            if(cmpSolution(newS, S)) S = newS;
            else if((double)(rand() % 10000) / 10000.0 < exp((double)(newS.obsSqrCnt - S.obsSqrCnt) / T)) S = newS;
        }
        T -= T_DEC;
    }

    return res;
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
            double eff = 100.0 * (double)res.obsSqrCnt / (double)upbOPT;
            /* Early termination
            if(eff + EPS > minEfficiency) tempCnt++; else tempCnt = 0;
            if(tempCnt >= 280) break;
            //*/
            minEfficiency = min(minEfficiency, eff);
            ///* Time calculation
            gettimeofday(&endTestTime, NULL);
            double testRunningTime = double(endTestTime.tv_sec - startTestTime.tv_sec) + double(endTestTime.tv_usec - startTestTime.tv_usec) * 1e-6;
            mapRunningTime += testRunningTime;
            //*/
            ///* Print test result
            //printSolution(res);
            cout << "T1 = " << originBmp.T1 << ", T2 = " << originBmp.T2 << ", C = " << originBmp.C << ", R = " << fixed << R << ", ";
            cout << "M = " << M << ", k = " << k << ": upbOPT = " << upbOPT << ", res = " << res.obsSqrCnt;
            cout << " | Eff: " << fixed << eff << ", minEff: " << fixed << minEfficiency;
            if(abso(eff - minEfficiency) <= EPS) cout << " (***)";
            cout << ", time: " << fixed << testRunningTime << " seconds";
            cout << endl;
            //*/
            /* Early termination
            if(upbOPT == originBmp.C - originBmp.unreachableSqrCnt) break;
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

iii getMandK(){
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 10 && abs(R - 0.500) < EPS) return iii(ii(2, 2), 6);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 10 && abs(R - 1.000) < EPS) return iii(ii(2, 1), 8);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 10 && abs(R - 2.000) < EPS) return iii(ii(1, 1), 8);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 20 && abs(R - 0.500) < EPS) return iii(ii(4, 1), 16);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 20 && abs(R - 1.000) < EPS) return iii(ii(2, 1), 13);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 20 && abs(R - 2.000) < EPS) return iii(ii(1, 1), 16);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 25 && abs(R - 0.500) < EPS) return iii(ii(4, 1), 19);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(3, 1), 17);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(2, 1), 22);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 50 && abs(R - 0.500) < EPS) return iii(ii(6, 1), 32);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 100 && abs(R - 1.000) < EPS) return iii(ii(11, 1), 61);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 100 && abs(R - 2.000) < EPS) return iii(ii(7, 1), 66);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 100 && abs(R - 3.000) < EPS) return iii(ii(3, 2), 72);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 10 && abs(R - 2.000) < EPS) return iii(ii(2, 2), 7);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 200 && abs(R - 1.000) < EPS) return iii(ii(5, 6), 124);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 200 && abs(R - 2.000) < EPS) return iii(ii(9, 1), 133);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 200 && abs(R - 3.000) < EPS) return iii(ii(3, 2), 138);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(4, 1), 13);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(4, 1), 16);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 25 && abs(R - 3.000) < EPS) return iii(ii(2, 2), 13);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 300 && abs(R - 1.000) < EPS) return iii(ii(6, 5), 194);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 300 && abs(R - 2.000) < EPS) return iii(ii(4, 4), 207);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 300 && abs(R - 3.000) < EPS) return iii(ii(5, 1), 200);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 400 && abs(R - 1.000) < EPS) return iii(ii(5, 9), 246);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 400 && abs(R - 2.000) < EPS) return iii(ii(4, 4), 273);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 50 && abs(R - 1.000) < EPS) return iii(ii(7, 1), 26);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 50 && abs(R - 2.000) < EPS) return iii(ii(4, 2), 35);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 50 && abs(R - 3.000) < EPS) return iii(ii(4, 1), 41);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 100 && abs(R - 1.000) < EPS) return iii(ii(5, 5), 58);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 100 && abs(R - 2.000) < EPS) return iii(ii(4, 3), 64);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 100 && abs(R - 3.000) < EPS) return iii(ii(3, 3), 74);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 200 && abs(R - 1.000) < EPS) return iii(ii(5, 7), 112);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 200 && abs(R - 2.000) < EPS) return iii(ii(3, 6), 120);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 200 && abs(R - 3.000) < EPS) return iii(ii(6, 1), 129);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(4, 2), 12);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(4, 1), 15);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 25 && abs(R - 3.000) < EPS) return iii(ii(3, 1), 16);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 400 && abs(R - 1.000) < EPS) return iii(ii(6, 8), 239);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 400 && abs(R - 2.000) < EPS) return iii(ii(4, 3), 234);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 400 && abs(R - 3.000) < EPS) return iii(ii(4, 2), 261);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 50 && abs(R - 1.000) < EPS) return iii(ii(9, 1), 30);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 50 && abs(R - 2.000) < EPS) return iii(ii(7, 1), 35);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 50 && abs(R - 3.000) < EPS) return iii(ii(3, 4), 36);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 100 && abs(R - 1.000) < EPS) return iii(ii(6, 5), 49);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 100 && abs(R - 2.000) < EPS) return iii(ii(5, 4), 62);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 100 && abs(R - 3.000) < EPS) return iii(ii(8, 1), 64);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 200 && abs(R - 1.000) < EPS) return iii(ii(9, 3), 106);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 200 && abs(R - 2.000) < EPS) return iii(ii(5, 6), 132);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 200 && abs(R - 3.000) < EPS) return iii(ii(9, 1), 122);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(4, 3), 14);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(7, 1), 17);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 25 && abs(R - 3.000) < EPS) return iii(ii(4, 1), 17);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 500 && abs(R - 1.000) < EPS) return iii(ii(8, 10), 282);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 500 && abs(R - 2.000) < EPS) return iii(ii(21, 1), 323);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 500 && abs(R - 3.000) < EPS) return iii(ii(11, 1), 303);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 50 && abs(R - 1.000) < EPS) return iii(ii(7, 2), 25);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 50 && abs(R - 2.000) < EPS) return iii(ii(5, 3), 33);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 50 && abs(R - 3.000) < EPS) return iii(ii(5, 1), 30);
}

void checkOneTest(){
    iii MandK = getMandK();
    M = MandK.first.first; k = MandK.first.second;
    int bestGreedyResult = MandK.second;
    struct timeval startTestTime, endTestTime;
    gettimeofday(&startTestTime, NULL);
    bmp = originBmp;
    int upbOPT = upperboundOPT();
    solution res = solveSOBP();
    //printSolution(res);
    double eff = 100.0 * (double)res.obsSqrCnt / (double)upbOPT;
    gettimeofday(&endTestTime, NULL);
    double testRunningTime = double(endTestTime.tv_sec - startTestTime.tv_sec) + double(endTestTime.tv_usec - startTestTime.tv_usec) * 1e-6;

    cout.precision(3);
    cout << "T1 = " << originBmp.T1 << ", T2 = " << originBmp.T2 << ", C = " << originBmp.C << ", R = " << fixed << R << ", ";
    cout << "M = " << M << ", k = " << k << ": upbOPT = " << upbOPT << ", res = " << res.obsSqrCnt << ", greedyRes = " << bestGreedyResult;
    cout << " | Eff: " << fixed << eff << ", greedyEff: " << fixed << 100.0 * (double)bestGreedyResult / (double)upbOPT << ", time: " << fixed << testRunningTime << " seconds" << endl;
}

void debug(){
    for(int i = 0; i < originBmp.N; i++){
        cout << "Bus #" << i + 1 << endl;
        for(int j = 0; j < (int)originBmp.buses[i].obsCritSqr.size(); j++){
            cout << "Critical point " << j + 1 << ":";
            for(int t = 0; t < (int)originBmp.buses[i].obsCritSqr[j].size(); t++){
                cout << "(" << originBmp.buses[i].obsCritSqr[j][t].first << ", " << originBmp.buses[i].obsCritSqr[j][t].second << "), ";
            }
            cout << endl;
        }
    }
}

int main(){
    srand(time(NULL));

    readInput();
    calculateUnreachableSqr();
    buildObsCritSqr();
    //debug();
    //fullCheck();
    checkOneTest();

    return 0;
}
