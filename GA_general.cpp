//Please comment all "Test old scenario" before running new tests
#include <bits/stdc++.h>
#include <sys/time.h>

using namespace std;

const double EPS = 1e-9;
const double INF = 1e9;
const double alpha = 1.0 - 1.0 / exp(1.0);

const int INITIAL_SOLUTIONS_FOR_EACH_BUS = 100;
const int G_BEST = 10;
const int G_AVG = 10;
const double P_MUTATION = 0.1;
int populationSize;
const bool USE_GREEDY_INITIAL_SOLUTION = false;

double sqr(double x){return x * x;}
double abso(double x){return (x < 0.0) ? -x : x;}

typedef pair<int, int> ii;
typedef pair<double, double> dd;
typedef pair<double, int> di;
typedef pair<dd, int> intv;
typedef pair<vector<int>, int> vii;
typedef pair<ii, int> iii;

struct busRoute{
    int busId;
    int numTurnPts[2]; //number of turning points
    vector<dd> turnPts[2]; //turning points
    int turnOnLimit[2];
    vector<double> critPtsLenFromStart[2]; // * length from start to a critical square
    vector<vector<ii> > obsCritSqr[2]; // * observable squares by a critical point
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
    vector<int> critPtsSet[2]; // * set of index of chosen critical points
};

struct solution{
    int obsSqrCnt; //observable critical squares
    vector<chosenBus> chosenBuses;
};

busMap bmp, originBmp;
double R;
int M, k;
bool fullCheckStarted = false;

vector<int> criticalSqrCoverCnt;

void readInput(){
    //freopen("Test/42x50_1000_0.50.txt", "r", stdin);

    cin >> bmp.T1 >> bmp.T2 >> R >> bmp.N;
    for(int i = 1; i <= bmp.N; i++){
        busRoute bus;
        bus.busId = i;
        cin >> bus.numTurnPts[0];
        for(int dir = 0; dir <= 1; dir++){
            cin >> bus.numTurnPts[dir];
            for(int j = 1; j <= bus.numTurnPts[dir]; j++){
                double x, y; cin >> x >> y;
                bus.turnPts[dir].push_back(dd(x, y));
            }
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

void calCoveredSqr(chosenBus* X){
    set<ii> obsSqr;
    for(int i = 0; i < (int)originBmp.buses.size(); i++){
        if(originBmp.buses[i].busId == (*X).busId){
            for(int dir = 0; dir <= 1; dir++){
                for(int j = 0; j < (int)(*X).critPtsSet[dir].size(); j++){
                    for(int t = 0; t < (int)originBmp.buses[i].obsCritSqr[dir][(*X).critPtsSet[dir][j]].size(); t++){
                        obsSqr.insert(originBmp.buses[i].obsCritSqr[dir][(*X).critPtsSet[dir][j]][t]);
                    }
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

chosenBus subOptimalSet(int busIdInCurMap){
    chosenBus cb;
    busRoute Y = bmp.buses[busIdInCurMap];
    cb.busId = Y.busId;

    //Calculate the observable intervals
    vector<int> reachableSqr[2];
    vector<double> criticalPts[2];
    for(int dir = 0; dir <= 1; dir++){
        for(int i = 0; i < bmp.C; i++){
            dd X = dd((double)bmp.criticalSqr[i].first, (double)bmp.criticalSqr[i].second);
            vector<double> L;
            double lenFromStart = 0.0;
            //If starting point of Y can observe X
            if(disPtsToSqr(Y.turnPts[dir][0], X) < R + EPS){
                L.push_back(0.0);
            }
            //Loop through all intervals in Y
            for(int j = 0; j < Y.numTurnPts[dir] - 1; j++){
                vector<double> tempIntersections = intsSegAndExtSqr(Y.turnPts[dir][j], Y.turnPts[dir][j + 1], X, lenFromStart);
                while(!tempIntersections.empty()){
                    L.push_back(tempIntersections.back());
                    tempIntersections.pop_back();
                }
                lenFromStart += disPtsToPts(Y.turnPts[dir][j], Y.turnPts[dir][j + 1]);
            }
            //If ending point of Y can observe X
            if(Y.numTurnPts[dir] > 1 && disPtsToSqr(Y.turnPts[dir][Y.numTurnPts[dir] - 1], X) < R + EPS){
                L.push_back(lenFromStart);
            }
            //Get observable interval
            if(!L.empty()){
                criticalSqrCoverCnt[i]++;
                reachableSqr[dir].push_back(i);
                sort(L.begin(), L.end());
                vector<double> tempL;
                tempL.push_back(L[0]);
                for(int j = 1; j < (int)L.size(); j++){
                    if(abso(L[j] - L[j - 1]) > EPS){
                        tempL.push_back(L[j]);
                    }
                }
                L = tempL;
                for(int j = 0; j < (int)L.size(); j++) criticalPts[dir].push_back(L[j]);
            }
        }
    }

    // If just checking for unreachable squares
    if(!fullCheckStarted) return cb;

    // Find sensor's turned on positions and observable squares
    vector<double> turnedOnPos[2];
    vector<bool> obsSqr(bmp.C, false);
    for(int dir = 0; dir <= 1; dir++) sort(criticalPts[dir].begin(), criticalPts[dir].end());
    vector<int> pickCntPerDir(2, 0);
    for(int i = 1; i <= Y.turnOnLimit[0] + Y.turnOnLimit[1]; i++){
        int pickDir = 0;
        double pickLen = -1;
        vector<int> maxCover;
        for(int dir = 0; dir <= 1; dir++){
            if(pickCntPerDir[dir] >= Y.turnOnLimit[dir]) continue;
            int lastTurnPts = 0;
            double lenFromStart = 0.0;
            for(int j = 0; j < (int)criticalPts[dir].size(); j++){
                while(lastTurnPts < Y.numTurnPts[dir] - 2
                   && criticalPts[dir][j] > lenFromStart + disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]) + EPS){
                    lenFromStart += disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]);
                    lastTurnPts++;
                }
                double len = criticalPts[dir][j] - lenFromStart;
                dd X;
                if(Y.numTurnPts[dir] == 1) X = dd((double)Y.turnPts[dir][0].first, (double)Y.turnPts[dir][0].second);
                else{
                    double dx = Y.turnPts[dir][lastTurnPts + 1].first - Y.turnPts[dir][lastTurnPts].first;
                    double dy = Y.turnPts[dir][lastTurnPts + 1].second - Y.turnPts[dir][lastTurnPts].second;
                    double sx = (dx + EPS < 0.0) ? -1.0 : 1.0;
                    double sy = (dy + EPS < 0.0) ? -1.0 : 1.0;
                    double a, b;
                    if(abso(dx) < EPS){
                        a = 0.0;
                        b = len;
                    }
                    else{
                        a = len / sqrt(1.0 + sqr(dy / dx));
                        b = sqrt(sqr(len) - sqr(a));
                    }
                    X.first = Y.turnPts[dir][lastTurnPts].first + sx * a;
                    X.second = Y.turnPts[dir][lastTurnPts].second + sy * b;
                }

                vector<int> tempCover;
                for(int t = 0; t < (int)reachableSqr[dir].size(); t++){
                    if(!obsSqr[reachableSqr[dir][t]]
                       && disPtsToSqr(X, dd((double)bmp.criticalSqr[reachableSqr[dir][t]].first, (double)bmp.criticalSqr[reachableSqr[dir][t]].second)) < R + EPS){
                        tempCover.push_back(reachableSqr[dir][t]);
                    }
                }
                if((int)tempCover.size() > (int)maxCover.size()){
                    maxCover = tempCover;
                    pickDir = dir;
                    pickLen = criticalPts[dir][j];
                }
            }
        }

        pickCntPerDir[pickDir]++;
        while(!maxCover.empty()){
            obsSqr[maxCover.back()] = true;
            maxCover.pop_back();
        }
        if(pickLen > -EPS) turnedOnPos[pickDir].push_back(pickLen);
    }
    for(int i = 0; i < bmp.C; i++) if(obsSqr[i]) cb.coveredSqr.push_back(bmp.criticalSqr[i]);

    // Trace critical points from the original bus map
    for(int dir = 0; dir <= 1; dir++){
        vector<bool> usedCritPts((int)bmp.buses[busIdInCurMap].critPtsLenFromStart[dir].size(), false);
        for(int i = 0; i < (int)turnedOnPos[dir].size(); i++){
            for(int j = 0; j < (int)bmp.buses[busIdInCurMap].critPtsLenFromStart[dir].size(); j++){
                if(usedCritPts[j]) continue;
                if(abs(turnedOnPos[dir][i] - bmp.buses[busIdInCurMap].critPtsLenFromStart[dir][j]) < EPS){
                    cb.critPtsSet[dir].push_back(j);
                    usedCritPts[j] = true;
                    break;
                }
            }
        }
    }
    calCoveredSqr(&cb);

    return cb;
}

bool subsetRelation(vector<int> X, vector<int> Y){
    for(int i = 0; i < (int)X.size(); i++){
        bool contained = false;
        for(int j = 0; j < (int)Y.size(); j++){
            if(X[i] == Y[j]){
                contained = true;
                break;
            }
        }
        if(!contained) return false;
    }
    return true;
}

int subOptimalSetSize(vector<vii> X, int turnOnLimit[2]){
    if(X.empty()) return 0;

    vector<int> pickCntPerDir(2, 0);
    vector<int> unionSet;
    for(int timer = 1; timer <= turnOnLimit[0] + turnOnLimit[1]; timer++){
        int maxNewElement = -1, chosenId = -1;
        for(int i = 0; i < (int)X.size(); i++){
            if(pickCntPerDir[X[i].second] >= turnOnLimit[X[i].second]) continue;
            int newElementCnt = 0;
            for(int j = 0; j < (int)X[i].first.size(); j++){
                bool exist = false;
                for(int t = 0; t < (int)unionSet.size(); t++){
                    if(X[i].first[j] == unionSet[t]){
                        exist = true;
                        break;
                    }
                }
                if(!exist) newElementCnt++;
            }
            if(newElementCnt > maxNewElement){
                chosenId = i;
                maxNewElement = newElementCnt;
            }
        }
        if(chosenId != -1){
            pickCntPerDir[X[chosenId].second]++;
            for(int i = 0; i < (int)X[chosenId].first.size(); i++){
                bool exist = false;
                for(int j = 0; j < (int)unionSet.size(); j++){
                    if(X[chosenId].first[i] == unionSet[j]){
                        exist = true;
                        break;
                    }
                }
                if(!exist) unionSet.push_back(X[chosenId].first[i]);
            }
        }
    }

    return (int)unionSet.size();
}

int upperboundOptimalSetSize(int busIdInCurMap){
    busRoute Y = bmp.buses[busIdInCurMap];

    //Calculate the observable intervals
    vector<int> reachableSqr[2];
    vector<double> criticalPts[2];
    for(int dir = 0; dir <= 1; dir++){
        for(int i = 0; i < bmp.C; i++){
            dd X = dd((double)bmp.criticalSqr[i].first, (double)bmp.criticalSqr[i].second);
            vector<double> L;
            double lenFromStart = 0.0;
            //If starting point of Y can observe X
            if(disPtsToSqr(Y.turnPts[dir][0], X) < R + EPS){
                L.push_back(0.0);
            }
            //Loop through all intervals in Y
            for(int j = 0; j < Y.numTurnPts[dir] - 1; j++){
                vector<double> tempIntersections = intsSegAndExtSqr(Y.turnPts[dir][j], Y.turnPts[dir][j + 1], X, lenFromStart);
                while(!tempIntersections.empty()){
                    L.push_back(tempIntersections.back());
                    tempIntersections.pop_back();
                }
                lenFromStart += disPtsToPts(Y.turnPts[dir][j], Y.turnPts[dir][j + 1]);
            }
            //If ending point of Y can observe X
            if(Y.numTurnPts[dir] > 1 && disPtsToSqr(Y.turnPts[dir][Y.numTurnPts[dir] - 1], X) < R + EPS){
                L.push_back(lenFromStart);
            }
            //Get observable interval
            if(!L.empty()){
                reachableSqr[dir].push_back(i);
                sort(L.begin(), L.end());
                vector<double> tempL;
                tempL.push_back(L[0]);
                for(int j = 1; j < (int)L.size(); j++){
                    if(abso(L[j] - L[j - 1]) > EPS){
                        tempL.push_back(L[j]);
                    }
                }
                L = tempL;
                for(int j = 0; j < (int)L.size(); j++) criticalPts[dir].push_back(L[j]);
            }
        }
    }

    // Find sensor's turned on positions and observable squares
    vector<vector<int> > coverSets;
    vector<vii> coverSetsWithDir;
    for(int dir = 0; dir <= 1; dir++){
        sort(criticalPts[dir].begin(), criticalPts[dir].end());
        int lastTurnPts = 0;
        double lenFromStart = 0.0;
        for(int i = 0; i < (int)criticalPts[dir].size(); i++){
            if(i > 0 && abso(criticalPts[dir][i] - criticalPts[dir][i-1]) < EPS) continue;
            while(lastTurnPts < Y.numTurnPts[dir] - 2
               && criticalPts[dir][i] > lenFromStart + disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]) + EPS){
                lenFromStart += disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]);
                lastTurnPts++;
            }
            double len = criticalPts[dir][i] - lenFromStart;
            dd X;
            if(Y.numTurnPts[dir] == 1) X = dd((double)Y.turnPts[dir][0].first, (double)Y.turnPts[dir][0].second);
            else{
                double dx = Y.turnPts[dir][lastTurnPts + 1].first - Y.turnPts[dir][lastTurnPts].first;
                double dy = Y.turnPts[dir][lastTurnPts + 1].second - Y.turnPts[dir][lastTurnPts].second;
                double sx = (dx + EPS < 0.0) ? -1.0 : 1.0;
                double sy = (dy + EPS < 0.0) ? -1.0 : 1.0;
                double a, b;
                if(abso(dx) < EPS){
                    a = 0.0;
                    b = len;
                }
                else{
                    a = len / sqrt(1.0 + sqr(dy / dx));
                    b = sqrt(sqr(len) - sqr(a));
                }
                X.first = Y.turnPts[dir][lastTurnPts].first + sx * a;
                X.second = Y.turnPts[dir][lastTurnPts].second + sy * b;
            }

            vector<int> tempCover;
            for(int t = 0; t < (int)reachableSqr[dir].size(); t++){
                if(disPtsToSqr(X, dd((double)bmp.criticalSqr[reachableSqr[dir][t]].first, (double)bmp.criticalSqr[reachableSqr[dir][t]].second)) < R + EPS){
                    tempCover.push_back(reachableSqr[dir][t]);
                }
            }
            coverSets.push_back(tempCover);
            coverSetsWithDir.push_back(vii(tempCover, dir));
        }
    }

    vector<int> coverCnt;
    for(int i = 0; i < (int)coverSets.size(); i++){
        bool ignoreSet = false;
        for(int j = 0; j < (int)coverSets.size(); j++){
            if(i == j) continue;
            if((int)coverSets[i].size() > (int)coverSets[j].size()) continue;
            if(subsetRelation(coverSets[i], coverSets[j])){
                if((int)coverSets[i].size() < (int)coverSets[j].size()
                   || (j < i && coverSets[i] == coverSets[j])){
                    ignoreSet = true;
                    break;
                }
            }
        }
        if(!ignoreSet) coverCnt.push_back((int)coverSets[i].size());
    }
    sort(coverCnt.begin(), coverCnt.end());

    int res = 0;
    for(int i = 1; i <= Y.turnOnLimit[0] + Y.turnOnLimit[1]; i++){
        if(coverCnt.empty()) break;
        res += coverCnt.back();
        coverCnt.pop_back();
    }

    return min(res, (int)ceil((double)subOptimalSetSize(coverSetsWithDir, Y.turnOnLimit) / alpha));
}

void buildObsCritSqrForBus(int busIdInCurMap){
    busRoute Y = originBmp.buses[busIdInCurMap];

    //Calculate the critical points on each path
    vector<double> criticalPts[2];
    for(int dir = 0; dir <= 1; dir++){
        for(int i = 0; i < originBmp.C; i++){
            dd X = dd((double)originBmp.criticalSqr[i].first, (double)originBmp.criticalSqr[i].second);
            vector<double> L;
            double lenFromStart = 0.0;
            //If starting point of Y can observe X
            if(disPtsToSqr(Y.turnPts[dir][0], X) < R + EPS){
                L.push_back(0.0);
            }
            //Loop through all intervals in Y
            for(int j = 0; j < Y.numTurnPts[dir] - 1; j++){
                vector<double> tempIntersections = intsSegAndExtSqr(Y.turnPts[dir][j], Y.turnPts[dir][j + 1], X, lenFromStart);
                while(!tempIntersections.empty()){
                    L.push_back(tempIntersections.back());
                    tempIntersections.pop_back();
                }
                lenFromStart += disPtsToPts(Y.turnPts[dir][j], Y.turnPts[dir][j + 1]);
            }
            //If ending point of Y can observe X
            if(Y.numTurnPts[dir] > 1 && disPtsToSqr(Y.turnPts[dir][Y.numTurnPts[dir] - 1], X) < R + EPS){
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
                for(int j = 0; j < (int)L.size(); j++) criticalPts[dir].push_back(L[j]);
            }
        }
    }

    // Calculate observable critical square set from each critical point
    for(int dir = 0; dir <= 1; dir++){
        sort(criticalPts[dir].begin(), criticalPts[dir].end());
        int lastTurnPts = 0;
        double lenFromStart = 0.0;
        for(int j = 0; j < (int)criticalPts[dir].size(); j++){
            while(lastTurnPts < Y.numTurnPts[dir] - 2
                && criticalPts[dir][j] > lenFromStart + disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]) + EPS){
                lenFromStart += disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]);
                lastTurnPts++;
            }
            double len = criticalPts[dir][j] - lenFromStart;
            dd X;
            if(Y.numTurnPts[dir] == 1) X = dd((double)Y.turnPts[dir][0].first, (double)Y.turnPts[dir][0].second);
            else{
                double dx = Y.turnPts[dir][lastTurnPts + 1].first - Y.turnPts[dir][lastTurnPts].first;
                double dy = Y.turnPts[dir][lastTurnPts + 1].second - Y.turnPts[dir][lastTurnPts].second;
                double sx = (dx + EPS < 0.0) ? -1.0 : 1.0;
                double sy = (dy + EPS < 0.0) ? -1.0 : 1.0;
                double a, b;
                if(abso(dx) < EPS){
                    a = 0.0;
                    b = len;
                }
                else{
                    a = len / sqrt(1.0 + sqr(dy / dx));
                    b = sqrt(sqr(len) - sqr(a));
                }
                X.first = Y.turnPts[dir][lastTurnPts].first + sx * a;
                X.second = Y.turnPts[dir][lastTurnPts].second + sy * b;
            }

            vector<ii> tempCover;
            for(int t = 0; t < (int)originBmp.criticalSqr.size(); t++){
                if(disPtsToSqr(X, dd((double)originBmp.criticalSqr[t].first, (double)originBmp.criticalSqr[t].second)) < R + EPS){
                    tempCover.push_back(originBmp.criticalSqr[t]);
                }
            }
            originBmp.buses[busIdInCurMap].critPtsLenFromStart[dir].push_back(criticalPts[dir][j]);
            originBmp.buses[busIdInCurMap].obsCritSqr[dir].push_back(tempCover);
        }
    }
}

void buildObsCritSqr(){
    for(int i = 0; i < (int)originBmp.buses.size(); i++){
        buildObsCritSqrForBus(i);
    }
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

    for(int i = 0; i < originBmp.N; i++) subOptimalSet(i);

    originBmp.unreachableSqrCnt = 0;
    for(int i = 0; i < originBmp.C; i++) originBmp.unreachableSqrCnt += (criticalSqrCoverCnt[i] == 0) ? 1 : 0;
    fullCheckStarted = true;
}

int upperboundOPT(){
    vector<int> cvSqr;
    for(int i = 0; i < bmp.N; i++) cvSqr.push_back(upperboundOptimalSetSize(i));
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
    for(int dir = 0; dir <= 1; dir++){
        vector<bool> pickedCritPts((int)originBmp.buses[index].obsCritSqr[dir].size(), false);
        for(int i = 1; i <= min(k, (int)originBmp.buses[index].obsCritSqr[dir].size()); i++){
            int cnt = rand() % ((int)originBmp.buses[index].obsCritSqr[dir].size() - i + 1) + 1;
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
        for(int i = 0; i < (int)pickedCritPts.size(); i++) if(pickedCritPts[i]) res.critPtsSet[dir].push_back(i);
    }
    calCoveredSqr(&res);
    return res;
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
            chosenBus Y = subOptimalSet(j);
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

vector<solution> buildPopulation(){
    populationSize = originBmp.N * INITIAL_SOLUTIONS_FOR_EACH_BUS;
    vector<solution> res;
    for(int i = 0; i < originBmp.N; i++){
        for(int timer = 1; timer <= INITIAL_SOLUTIONS_FOR_EACH_BUS; timer++){
            solution S;
            S.obsSqrCnt = 0;
            S.chosenBuses.push_back(createRandomizedChosenBus(originBmp.buses[i].busId));
            vector<bool> pickedBuses(originBmp.N, false);
            pickedBuses[i] = true;
            for(int t = 1; t <= M - 1; t++){
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
            for(int j = 0; j < (int)pickedBuses.size(); j++) if(pickedBuses[j] && j != i) S.chosenBuses.push_back(createRandomizedChosenBus(originBmp.buses[j].busId));
            calSolutionStats(&S);
            res.push_back(S);
        }
    }
    return res;
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

double calAvgObsSqrCnt(vector<solution> population){
    int sum = 0;
    for(int i = 0; i < populationSize; i++) sum += population[i].obsSqrCnt;
    return (double)sum / (double)populationSize;
}

solution crossover(solution A, solution B){
    solution res;
    vector<int> busIndexesA(originBmp.N + 1, -1);
    vector<int> busIndexesB(originBmp.N + 1, -1);
    for(int i = 0; i < (int)A.chosenBuses.size(); i++) busIndexesA[A.chosenBuses[i].busId] = i;
    for(int i = 0; i < (int)B.chosenBuses.size(); i++) busIndexesB[B.chosenBuses[i].busId] = i;
    vector<chosenBus> freeToChooseBuses;
    for(int i = 1; i <= originBmp.N; i++){
        if(busIndexesA[i] != -1 && busIndexesB[i] != -1){
            int side = rand() % 2;
            if(side == 0) res.chosenBuses.push_back(A.chosenBuses[busIndexesA[i]]);
            else res.chosenBuses.push_back(B.chosenBuses[busIndexesB[i]]);
        }
        else if(busIndexesA[i] != -1) freeToChooseBuses.push_back(A.chosenBuses[busIndexesA[i]]);
        else if(busIndexesB[i] != -1) freeToChooseBuses.push_back(B.chosenBuses[busIndexesB[i]]);
    }
    int remainBusCnt = (int)freeToChooseBuses.size() / 2;
    while(remainBusCnt > 0){
        int id = rand() % (int)freeToChooseBuses.size();
        res.chosenBuses.push_back(freeToChooseBuses[id]);
        freeToChooseBuses.erase(freeToChooseBuses.begin() + id);
        remainBusCnt--;
    }
    calSolutionStats(&res);
    return res;
}

vector<solution> individualSelection(vector<solution> population){
    vector<solution> newPopulation;
    for(int timer = 1; timer <= populationSize; timer++){
        ///* This line is for "greedy individual selection", comment it for "probabilistic individual selection"
        newPopulation.push_back(population[timer - 1]); continue;
        //*/
        int sum = 0;
        for(int i = 0; i < (int)population.size(); i++) sum += population[i].obsSqrCnt;
        int id;
        if(sum == 0) id = rand() % (int)population.size();
        else{
            int cnt = rand() % sum + 1;
            sum = 0; id = -1;
            while(sum < cnt){
                id++;
                sum += population[id].obsSqrCnt;
            }
        }
        newPopulation.push_back(population[id]);
        population.erase(population.begin() + id);
    }
    sort(newPopulation.begin(), newPopulation.end(), cmpSolution);
    return newPopulation;
}

solution solveSOBP(){
    /// GA all generalAndSpecial algorithm
    vector<solution> population = buildPopulation();
    sort(population.begin(), population.end(), cmpSolution);
    if(USE_GREEDY_INITIAL_SOLUTION) population[(int)population.size() - 1] = createGreedyInitialSolution();
    sort(population.begin(), population.end(), cmpSolution);
    solution res = population[0];
    int bestObsSqrCntIdle = 0, avgObsSqrCntIdle = 0;
    double bestAvgObsSqrCnt = calAvgObsSqrCnt(population);

    //int generationCnt = 0;
    //cout << "Generation " << generationCnt << " (POPULATION_SIZE = " << populationSize << "): bestRes = " << res.obsSqrCnt << ", bestAvg = " << bestAvgObsSqrCnt << endl;

    while(bestObsSqrCntIdle < G_BEST || avgObsSqrCntIdle < G_AVG){
        vector<solution> newPopulation = population;

        vector<int> sumObsSqrCnt(populationSize, 0);
        sumObsSqrCnt[0] = population[0].obsSqrCnt;
        for(int i = 1; i < populationSize; i++) sumObsSqrCnt[i] = sumObsSqrCnt[i - 1] + population[i].obsSqrCnt;
        for(int i = 0; i < populationSize; i++){
            /// Crossover
            if(sumObsSqrCnt[populationSize - 1] > 0){
                int cnt = rand() % sumObsSqrCnt[populationSize - 1] + 1;
                int le = 0, ri = populationSize - 1;
                while(le < ri){
                    int mid = (le + ri) / 2;
                    if(sumObsSqrCnt[mid] < cnt) le = mid + 1;
                    else ri = mid;
                }
                if(le != i){
                    newPopulation.push_back(crossover(population[i], population[le]));
                }
                else if(i > 0){
                    newPopulation.push_back(crossover(population[i], population[i - 1]));
                }
                else if(i < populationSize - 1){
                    newPopulation.push_back(crossover(population[i], population[i + 1]));
                }
            }
            else{
                newPopulation.push_back(crossover(population[i], population[rand() % populationSize]));
            }

            /// Mutation
            if((double)(rand() % 100) / 100.0 < P_MUTATION){
                // Swap a bus
                vector<bool> usedBus(originBmp.N + 1, false);
                for(int j = 0; j < M; j++) usedBus[population[i].chosenBuses[j].busId] = true;
                int busToRemove = rand() % M;
                usedBus[population[i].chosenBuses[busToRemove].busId] = false;
                solution newSolution = population[i];
                newSolution.chosenBuses.erase(newSolution.chosenBuses.begin() + busToRemove);
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
                newSolution.chosenBuses.push_back(createRandomizedChosenBus(busToAdd));
                calSolutionStats(&newSolution);
                if(population[i].chosenBuses[busToRemove].busId != busToAdd) newPopulation.push_back(newSolution);

                // Swap a critical point
                int busToSwap = rand() % M, busToSwapIndex = -1;
                for(int j = 0; j < originBmp.N; j++){
                    if(originBmp.buses[j].busId == population[i].chosenBuses[busToSwap].busId){
                        busToSwapIndex = j;
                        break;
                    }
                }
                int dirToSwap = rand() % 2;
                if((int)originBmp.buses[busToSwapIndex].obsCritSqr[dirToSwap].size() > 0){
                    vector<bool> usedCritPts((int)originBmp.buses[busToSwapIndex].obsCritSqr[dirToSwap].size(), false);
                    for(int j = 0; j < (int)population[i].chosenBuses[busToSwap].critPtsSet[dirToSwap].size(); j++) usedCritPts[population[i].chosenBuses[busToSwap].critPtsSet[dirToSwap][j]] = true;
                    if((int)population[i].chosenBuses[busToSwap].critPtsSet[dirToSwap].size() > 0){
                        int critPtsToRemove = rand() % (int)population[i].chosenBuses[busToSwap].critPtsSet[dirToSwap].size();
                        usedCritPts[population[i].chosenBuses[busToSwap].critPtsSet[dirToSwap][critPtsToRemove]] = false;
                        newSolution = population[i];
                        newSolution.chosenBuses[busToSwap].critPtsSet[dirToSwap].erase(newSolution.chosenBuses[busToSwap].critPtsSet[dirToSwap].begin() + critPtsToRemove);
                        int cntCritPtsToAdd = rand() % ((int)originBmp.buses[busToSwapIndex].obsCritSqr[dirToSwap].size() - (int)newSolution.chosenBuses[busToSwap].critPtsSet[dirToSwap].size()) + 1, critPtsToAdd = -1;
                        for(int j = 0; j < (int)originBmp.buses[busToSwapIndex].obsCritSqr[dirToSwap].size(); j++){
                            if(!usedCritPts[j]){
                                cntCritPtsToAdd--;
                                if(cntCritPtsToAdd == 0){
                                    critPtsToAdd = j;
                                    break;
                                }
                            }
                        }
                        newSolution.chosenBuses[busToSwap].critPtsSet[dirToSwap].push_back(critPtsToAdd);
                        calCoveredSqr(&(newSolution.chosenBuses[busToSwap]));
                        calSolutionStats(&newSolution);
                        if(critPtsToAdd != population[i].chosenBuses[busToSwap].critPtsSet[dirToSwap][critPtsToRemove]) newPopulation.push_back(newSolution);
                    }
                }
            }
        }

        // Update terminate conditions
        //generationCnt++;
        //cout << "Generation " << generationCnt << " (candidates = " << (int)newPopulation.size();

        sort(newPopulation.begin(), newPopulation.end(), cmpSolution);
        if(cmpSolution(newPopulation[0], res)){
            res = newPopulation[0];
            bestObsSqrCntIdle = 0;
        }
        else bestObsSqrCntIdle++;

        newPopulation = individualSelection(newPopulation);
        double curAvgObsSqrCnt = calAvgObsSqrCnt(newPopulation);
        if(curAvgObsSqrCnt > bestAvgObsSqrCnt + EPS){
            bestAvgObsSqrCnt = curAvgObsSqrCnt;
            avgObsSqrCntIdle = 0;
        }
        else avgObsSqrCntIdle++;

        // Set new population as current population
        population = newPopulation;

        //cout << "): bestRes = " << res.obsSqrCnt << ", bestAvg = " << bestAvgObsSqrCnt << endl;
        //printSolution(population[0]);
        //printSolution(population[populationSize - 1]);
    }
    
    return res;
}

vector<int> createDeltaSet(){
    vector<int> res;
    res.push_back(0);
    /* Additional values for delta
    res.push_back(5);
    res.push_back(10);
    res.push_back(20);
    res.push_back(50);
    res.push_back(originBmp.C/16);
    res.push_back(originBmp.C/8);
    res.push_back(originBmp.C/4);
    res.push_back(originBmp.C/2);
    res.push_back(originBmp.C);
    //*/
    return res;
}

void createTurnOnLimit(int delta){
    int leftLimit = max(0, k - delta);
    int rightLimit = min(originBmp.C, k + delta);
    for(int i = 0; i < originBmp.N; i++){
        for(int j = 0; j <= 1; j++){
            originBmp.buses[i].turnOnLimit[j] = leftLimit + rand() % (rightLimit - leftLimit + 1);
        }
    }
}

// Check all cases of M and k
void fullCheck(){
    fullCheckStarted = true;
    cout.precision(3);
    ///* Time calculation
    double mapRunningTime = 0.0;
    //*/
    /* Early termination *
    int tempCnt = 0;
    //*/
    int test_cnt = 0;
    double minEfficiency = 100.0;
    for(k = 1; k <= originBmp.C; k++){
        vector<int> deltaSet = createDeltaSet();
        for(int i = 0; i < (int)deltaSet.size(); i++){
            int delta = deltaSet[i];
            createTurnOnLimit(delta);
            /* Early termination
            int last_res = 0;
            int last_upbOPT = 0;
            //*/
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
                /* Early termination *
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
                cout << "M = " << M << ", k = " << k << ", delta = " << delta << ": upbOPT = " << upbOPT << ", res = " << res.obsSqrCnt;
                cout << " | Eff: " << fixed << eff << ", minEff: " << fixed << minEfficiency;
                if(abso(eff - minEfficiency) <= EPS) cout << " (***)";
                cout << ", time: " << fixed << testRunningTime << " seconds";
                cout << endl;
                //*/

                /* Early termination
                if(last_res == res.obsSqrCnt && last_upbOPT == upbOPT) break;
                last_res = res.obsSqrCnt;
                last_upbOPT = upbOPT;
                //*/

                /* Early termination
                if(upbOPT == originBmp.C - originBmp.unreachableSqrCnt) break;
                //*/
            }
            /* Early termination *
            if(tempCnt >= 280) break;
            //*/
        }
        /* Early termination *
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
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 10 && abs(R - 0.500) < EPS) return iii(ii(1, 2), 7);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 10 && abs(R - 1.000) < EPS) return iii(ii(1, 2), 6);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 10 && abs(R - 2.000) < EPS) return iii(ii(1, 2), 8);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 20 && abs(R - 0.500) < EPS) return iii(ii(2, 2), 12);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 20 && abs(R - 1.000) < EPS) return iii(ii(1, 3), 11);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 20 && abs(R - 2.000) < EPS) return iii(ii(1, 2), 15);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 25 && abs(R - 0.500) < EPS) return iii(ii(1, 5), 12);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(2, 1), 15);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(1, 2), 18);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 30 && abs(R - 0.500) < EPS) return iii(ii(2, 4), 19);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 30 && abs(R - 1.000) < EPS) return iii(ii(1, 4), 17);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 30 && abs(R - 2.000) < EPS) return iii(ii(1, 2), 21);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 40 && abs(R - 0.500) < EPS) return iii(ii(2, 3), 24);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 40 && abs(R - 1.000) < EPS) return iii(ii(1, 4), 20);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 40 && abs(R - 2.000) < EPS) return iii(ii(1, 2), 28);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 50 && abs(R - 0.500) < EPS) return iii(ii(1, 5), 19);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 50 && abs(R - 1.000) < EPS) return iii(ii(1, 5), 24);
    if(originBmp.T1 == 10 && originBmp.T2 == 12 && originBmp.C == 50 && abs(R - 2.000) < EPS) return iii(ii(1, 2), 37);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 100 && abs(R - 0.500) < EPS) return iii(ii(4, 7), 52);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 100 && abs(R - 1.000) < EPS) return iii(ii(2, 9), 39);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 100 && abs(R - 2.000) < EPS) return iii(ii(2, 4), 57);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 100 && abs(R - 3.000) < EPS) return iii(ii(1, 3), 40);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 10 && abs(R - 0.500) < EPS) return iii(ii(3, 2), 7);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 10 && abs(R - 1.000) < EPS) return iii(ii(3, 1), 6);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 10 && abs(R - 2.000) < EPS) return iii(ii(1, 1), 3);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 10 && abs(R - 3.000) < EPS) return iii(ii(2, 1), 8);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 200 && abs(R - 0.500) < EPS) return iii(ii(4, 12), 104);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 200 && abs(R - 1.000) < EPS) return iii(ii(2, 12), 77);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 200 && abs(R - 2.000) < EPS) return iii(ii(2, 10), 105);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 200 && abs(R - 3.000) < EPS) return iii(ii(1, 5), 95);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 25 && abs(R - 0.500) < EPS) return iii(ii(6, 1), 15);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(5, 1), 15);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(3, 1), 14);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 25 && abs(R - 3.000) < EPS) return iii(ii(1, 5), 15);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 300 && abs(R - 0.500) < EPS) return iii(ii(3, 20), 114);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 300 && abs(R - 1.000) < EPS) return iii(ii(2, 12), 120);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 300 && abs(R - 2.000) < EPS) return iii(ii(2, 6), 159);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 300 && abs(R - 3.000) < EPS) return iii(ii(1, 9), 122);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 400 && abs(R - 0.500) < EPS) return iii(ii(3, 19), 148);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 400 && abs(R - 1.000) < EPS) return iii(ii(2, 16), 157);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 400 && abs(R - 2.000) < EPS) return iii(ii(2, 6), 210);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 400 && abs(R - 3.000) < EPS) return iii(ii(1, 6), 181);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 50 && abs(R - 0.500) < EPS) return iii(ii(3, 6), 23);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 50 && abs(R - 1.000) < EPS) return iii(ii(3, 5), 32);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 50 && abs(R - 2.000) < EPS) return iii(ii(2, 3), 31);
    if(originBmp.T1 == 25 && originBmp.T2 == 30 && originBmp.C == 50 && abs(R - 3.000) < EPS) return iii(ii(1, 5), 26);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 100 && abs(R - 0.500) < EPS) return iii(ii(5, 5), 50);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 100 && abs(R - 1.000) < EPS) return iii(ii(10, 1), 49);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 100 && abs(R - 2.000) < EPS) return iii(ii(2, 4), 50);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 100 && abs(R - 3.000) < EPS) return iii(ii(1, 6), 42);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 10 && abs(R - 0.500) < EPS) return iii(ii(3, 1), 7);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 10 && abs(R - 1.000) < EPS) return iii(ii(3, 2), 7);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 10 && abs(R - 2.000) < EPS) return iii(ii(1, 2), 6);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 10 && abs(R - 3.000) < EPS) return iii(ii(1, 3), 5);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 200 && abs(R - 0.500) < EPS) return iii(ii(4, 11), 93);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 200 && abs(R - 1.000) < EPS) return iii(ii(3, 11), 91);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 200 && abs(R - 2.000) < EPS) return iii(ii(2, 10), 100);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 200 && abs(R - 3.000) < EPS) return iii(ii(1, 6), 78);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 25 && abs(R - 0.500) < EPS) return iii(ii(5, 2), 16);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(5, 1), 12);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(3, 2), 19);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 25 && abs(R - 3.000) < EPS) return iii(ii(3, 1), 17);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 300 && abs(R - 0.500) < EPS) return iii(ii(4, 14), 127);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 300 && abs(R - 1.000) < EPS) return iii(ii(3, 12), 123);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 300 && abs(R - 2.000) < EPS) return iii(ii(2, 8), 154);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 300 && abs(R - 3.000) < EPS) return iii(ii(2, 5), 178);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 400 && abs(R - 0.500) < EPS) return iii(ii(4, 17), 159);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 400 && abs(R - 1.000) < EPS) return iii(ii(3, 14), 175);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 400 && abs(R - 2.000) < EPS) return iii(ii(2, 11), 187);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 400 && abs(R - 3.000) < EPS) return iii(ii(1, 9), 153);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 500 && abs(R - 0.500) < EPS) return iii(ii(3, 20), 172);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 500 && abs(R - 1.000) < EPS) return iii(ii(3, 12), 231);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 500 && abs(R - 2.000) < EPS) return iii(ii(2, 9), 241);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 500 && abs(R - 3.000) < EPS) return iii(ii(2, 5), 284);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 50 && abs(R - 0.500) < EPS) return iii(ii(4, 6), 26);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 50 && abs(R - 1.000) < EPS) return iii(ii(2, 6), 22);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 50 && abs(R - 2.000) < EPS) return iii(ii(2, 3), 24);
    if(originBmp.T1 == 30 && originBmp.T2 == 36 && originBmp.C == 50 && abs(R - 3.000) < EPS) return iii(ii(1, 4), 23);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 1000 && abs(R - 0.500) < EPS) return iii(ii(4, 28), 299);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 1000 && abs(R - 1.000) < EPS) return iii(ii(3, 26), 348);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 1000 && abs(R - 2.000) < EPS) return iii(ii(3, 11), 494);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 1000 && abs(R - 3.000) < EPS) return iii(ii(2, 20), 495);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 100 && abs(R - 0.500) < EPS) return iii(ii(6, 6), 45);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 100 && abs(R - 1.000) < EPS) return iii(ii(5, 7), 52);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 100 && abs(R - 2.000) < EPS) return iii(ii(3, 6), 55);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 100 && abs(R - 3.000) < EPS) return iii(ii(2, 8), 49);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 10 && abs(R - 0.500) < EPS) return iii(ii(1, 1), 2);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 10 && abs(R - 1.000) < EPS) return iii(ii(4, 1), 8);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 10 && abs(R - 2.000) < EPS) return iii(ii(2, 2), 7);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 10 && abs(R - 3.000) < EPS) return iii(ii(2, 1), 5);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 200 && abs(R - 0.500) < EPS) return iii(ii(8, 5), 96);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 200 && abs(R - 1.000) < EPS) return iii(ii(4, 11), 92);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 200 && abs(R - 2.000) < EPS) return iii(ii(2, 10), 79);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 200 && abs(R - 3.000) < EPS) return iii(ii(2, 6), 97);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 25 && abs(R - 0.500) < EPS) return iii(ii(2, 1), 5);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 25 && abs(R - 1.000) < EPS) return iii(ii(4, 3), 14);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 25 && abs(R - 2.000) < EPS) return iii(ii(6, 1), 19);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 25 && abs(R - 3.000) < EPS) return iii(ii(1, 3), 7);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 300 && abs(R - 0.500) < EPS) return iii(ii(6, 10), 114);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 300 && abs(R - 1.000) < EPS) return iii(ii(3, 14), 104);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 300 && abs(R - 2.000) < EPS) return iii(ii(3, 9), 150);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 300 && abs(R - 3.000) < EPS) return iii(ii(2, 9), 142);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 400 && abs(R - 0.500) < EPS) return iii(ii(5, 16), 144);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 400 && abs(R - 1.000) < EPS) return iii(ii(3, 17), 147);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 400 && abs(R - 2.000) < EPS) return iii(ii(2, 14), 156);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 400 && abs(R - 3.000) < EPS) return iii(ii(2, 7), 187);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 500 && abs(R - 0.500) < EPS) return iii(ii(6, 17), 209);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 500 && abs(R - 1.000) < EPS) return iii(ii(4, 16), 210);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 500 && abs(R - 2.000) < EPS) return iii(ii(2, 15), 193);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 500 && abs(R - 3.000) < EPS) return iii(ii(2, 10), 237);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 50 && abs(R - 0.500) < EPS) return iii(ii(5, 3), 21);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 50 && abs(R - 1.000) < EPS) return iii(ii(6, 3), 29);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 50 && abs(R - 2.000) < EPS) return iii(ii(3, 4), 25);
    if(originBmp.T1 == 42 && originBmp.T2 == 50 && originBmp.C == 50 && abs(R - 3.000) < EPS) return iii(ii(4, 1), 26);
}

void checkOneTest(){
    iii MandK = getMandK();
    M = MandK.first.first; k = MandK.first.second;
    int bestGreedyResult = MandK.second;
    struct timeval startTestTime, endTestTime;
    gettimeofday(&startTestTime, NULL);
    createTurnOnLimit(0);
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
    for(int i = 0; i < (int)originBmp.buses.size(); i++){
        cout << "Bus #" << i + 1 << endl;
        for(int dir = 0; dir <= 1; dir++){
            cout << "Direction " << dir + 1 << endl;
            for(int j = 0; j < (int)originBmp.buses[i].obsCritSqr[dir].size(); j++){
                cout << "Critical point " << j + 1 << ":";
                for(int t = 0; t < (int)originBmp.buses[i].obsCritSqr[dir][j].size(); t++){
                    cout << "(" << originBmp.buses[i].obsCritSqr[dir][j][t].first << ", " << originBmp.buses[i].obsCritSqr[dir][j][t].second << "), ";
                }
                cout << endl;
            }
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
