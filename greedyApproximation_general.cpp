//Please comment all "Test old scenario" before running new tests
#include <bits/stdc++.h>
#include <sys/time.h>

using namespace std;

const double eps = 1e-9;
const double INF = 1e9;
const double alpha = 1.0 - 1.0 / exp(1.0);

double sqr(double x){return x * x;}
double abso(double x){return (x < 0.0) ? -x : x;}

typedef pair<int, int> ii;
typedef pair<double, double> dd;
typedef pair<double, int> di;
typedef pair<dd, int> intv;
typedef pair<vector<int>, int> vii;

struct busRoute{
    int busId;
    int numTurnPts[2]; //number of turning points
    vector<dd> turnPts[2]; //turning points
    int turnOnLimit[2];
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
            if(disPtsToSqr(Y.turnPts[dir][0], X) < R + eps){
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
            if(Y.numTurnPts[dir] > 1 && disPtsToSqr(Y.turnPts[dir][Y.numTurnPts[dir] - 1], X) < R + eps){
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
                    if(abso(L[j] - L[j - 1]) > eps){
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
    vector<double> turnedOnPos;
    vector<bool> obsSqr(bmp.C, false);
    for(int dir = 0; dir <= 1; dir++) sort(criticalPts[dir].begin(), criticalPts[dir].end());
    vector<int> pickCntPerDir(2, 0);
    for(int i = 1; i <= Y.turnOnLimit[0] + Y.turnOnLimit[1]; i++){
        int pickDir = 0;
        vector<int> maxCover;
        for(int dir = 0; dir <= 1; dir++){
            if(pickCntPerDir[dir] >= Y.turnOnLimit[dir]) continue;
            int lastTurnPts = 0;
            double lenFromStart = 0.0;
            for(int j = 0; j < (int)criticalPts[dir].size(); j++){
                while(lastTurnPts < Y.numTurnPts[dir] - 2
                   && criticalPts[dir][j] > lenFromStart + disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]) + eps){
                    lenFromStart += disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]);
                    lastTurnPts++;
                }
                double len = criticalPts[dir][j] - lenFromStart;
                dd X;
                if(Y.numTurnPts[dir] == 1) X = dd((double)Y.turnPts[dir][0].first, (double)Y.turnPts[dir][0].second);
                else{
                    double dx = Y.turnPts[dir][lastTurnPts + 1].first - Y.turnPts[dir][lastTurnPts].first;
                    double dy = Y.turnPts[dir][lastTurnPts + 1].second - Y.turnPts[dir][lastTurnPts].second;
                    double sx = (dx + eps < 0.0) ? -1.0 : 1.0;
                    double sy = (dy + eps < 0.0) ? -1.0 : 1.0;
                    double a, b;
                    if(abso(dx) < eps){
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
                       && disPtsToSqr(X, dd((double)bmp.criticalSqr[reachableSqr[dir][t]].first, (double)bmp.criticalSqr[reachableSqr[dir][t]].second)) < R + eps){
                        tempCover.push_back(reachableSqr[dir][t]);
                    }
                }
                if((int)tempCover.size() > (int)maxCover.size()){
                    maxCover = tempCover;
                    pickDir = dir;
                }
            }
        }

        pickCntPerDir[pickDir]++;
        while(!maxCover.empty()){
            obsSqr[maxCover.back()] = true;
            maxCover.pop_back();
        }
    }
    for(int i = 0; i < bmp.C; i++) if(obsSqr[i]) cb.coveredSqr.push_back(bmp.criticalSqr[i]);

    return cb;
}

bool checkDuplicatePaths(int busId){
    if(bmp.buses[busId].numTurnPts[0] != bmp.buses[busId].numTurnPts[1]) return false;
    for(int i = 0; i < bmp.buses[busId].numTurnPts[0] / 2; i++){
        if(abso(bmp.buses[busId].turnPts[0][i].first - bmp.buses[busId].turnPts[1][bmp.buses[busId].numTurnPts[1] - 1 - i].first) > eps
           || abso(bmp.buses[busId].turnPts[0][i].second - bmp.buses[busId].turnPts[1][bmp.buses[busId].numTurnPts[1] - 1 - i].second) > eps){
            return false;
        }
    }
    return true;
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
            if(disPtsToSqr(Y.turnPts[dir][0], X) < R + eps){
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
            if(Y.numTurnPts[dir] > 1 && disPtsToSqr(Y.turnPts[dir][Y.numTurnPts[dir] - 1], X) < R + eps){
                L.push_back(lenFromStart);
            }
            //Get observable interval
            if(!L.empty()){
                reachableSqr[dir].push_back(i);
                sort(L.begin(), L.end());
                vector<double> tempL;
                tempL.push_back(L[0]);
                for(int j = 1; j < (int)L.size(); j++){
                    if(abso(L[j] - L[j - 1]) > eps){
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
            if(i > 0 && abso(criticalPts[dir][i] - criticalPts[dir][i-1]) < eps) continue;
            while(lastTurnPts < Y.numTurnPts[dir] - 2
               && criticalPts[dir][i] > lenFromStart + disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]) + eps){
                lenFromStart += disPtsToPts(Y.turnPts[dir][lastTurnPts], Y.turnPts[dir][lastTurnPts + 1]);
                lastTurnPts++;
            }
            double len = criticalPts[dir][i] - lenFromStart;
            dd X;
            if(Y.numTurnPts[dir] == 1) X = dd((double)Y.turnPts[dir][0].first, (double)Y.turnPts[dir][0].second);
            else{
                double dx = Y.turnPts[dir][lastTurnPts + 1].first - Y.turnPts[dir][lastTurnPts].first;
                double dy = Y.turnPts[dir][lastTurnPts + 1].second - Y.turnPts[dir][lastTurnPts].second;
                double sx = (dx + eps < 0.0) ? -1.0 : 1.0;
                double sy = (dy + eps < 0.0) ? -1.0 : 1.0;
                double a, b;
                if(abso(dx) < eps){
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
                if(disPtsToSqr(X, dd((double)bmp.criticalSqr[reachableSqr[dir][t]].first, (double)bmp.criticalSqr[reachableSqr[dir][t]].second)) < R + eps){
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

vector<ii> fullSet(int busIdInCurMap){
    vector<ii> res;
    busRoute Y = bmp.buses[busIdInCurMap];

    //Calculate the observable intervals
    for(int i = 0; i < bmp.C; i++){
        dd X = dd((double)bmp.criticalSqr[i].first, (double)bmp.criticalSqr[i].second);
        for(int dir = 0; dir <= 1; dir++){
            vector<double> L;
            double lenFromStart = 0.0;
            //If starting point of Y can observe X
            if(disPtsToSqr(Y.turnPts[dir][0], X) <= R + eps){
                L.push_back(0.0);
            }
            //Loop through all intervals in Y
            for(int j = 0; j < Y.numTurnPts[dir] - 1; j++){
                vector<double> tempIntersections = intsSegAndExtSqr(Y.turnPts[dir][j], Y.turnPts[dir][j+1], X, lenFromStart);
                while(!tempIntersections.empty()){
                    L.push_back(tempIntersections.back());
                    tempIntersections.pop_back();
                }
                lenFromStart += disPtsToPts(Y.turnPts[dir][j], Y.turnPts[dir][j+1]);
            }
            //If ending point of Y can observe X
            if(Y.numTurnPts[dir] > 1 && disPtsToSqr(Y.turnPts[dir][Y.numTurnPts[dir] - 1], X) <= R + eps){
                L.push_back(lenFromStart);
            }
            //Get observable interval
            if(!L.empty()){
                res.push_back(bmp.criticalSqr[i]);
                break;
            }
        }
    }

    return res;
}

void calculateUnreachableSqr(){
    fill(criticalSqrCoverCnt.begin(), criticalSqrCoverCnt.end(), 0);

    for(int i = 0; i < bmp.N; i++) subOptimalSet(i);

    originBmp.unreachableSqr = 0;
    for(int i = 0; i < bmp.C; i++) originBmp.unreachableSqr += (criticalSqrCoverCnt[i] == 0) ? 1 : 0;
}

int upperboundOPT(){
    vector<int> cvSqr;
    for(int i = 0; i < bmp.N; i++) cvSqr.push_back(upperboundOptimalSetSize(i));
    sort(cvSqr.begin(), cvSqr.end());
    int upbOPT = 0;
    for(int i = bmp.N - 1; i >= bmp.N - M; i--){
        upbOPT += cvSqr[i];
    }
    return min(upbOPT, bmp.C - bmp.unreachableSqr);
}

int maxCoveredSqr(){
    fill(criticalSqrCoverCnt.begin(), criticalSqrCoverCnt.end(), 0);

    for(int i = 0; i < bmp.N; i++) subOptimalSet(i);

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
            chosenBus Y = subOptimalSet(j);
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
            ///* Early termination
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
                double eff = 100.0 * (double)res.obsSqr / (double)upbOPT;
                /* Early termination *
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
                cout << "M = " << M << ", k = " << k << ", delta = " << delta << ": upbOPT = " << upbOPT << ", res = " << res.obsSqr;
                cout << " | Eff: " << fixed << eff << ", minEff: " << fixed << minEfficiency;
                if(abso(eff - minEfficiency) <= eps) cout << " (***)";
                cout << ", time: " << fixed << testRunningTime << " seconds";
                cout << endl;
                //*/

                ///* Early termination
                if(last_res == res.obsSqr && last_upbOPT == upbOPT) break;
                last_res = res.obsSqr;
                last_upbOPT = upbOPT;
                //*/

                /* Early termination
                if(upbOPT == originBmp.C - originBmp.unreachableSqr) break;
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

int main(){
    srand(time(NULL));
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
