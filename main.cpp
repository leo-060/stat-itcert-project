#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;
vector<vector<double>> identity_matrix(int N){
    vector temp (N, vector<double>(N, 0));
    for (int i=0;i<N;++i)
        temp[i][i] = 1;
    return temp;
}
vector<double> operator+(const vector<double> &a, const vector<double> &b){ // row addition
    vector<double> temp;
    for (int i=0;i<a.size();++i)
        temp.push_back(a[i] + b[i]);
    return temp;
}
vector<double> operator*(const double k, const vector<double> &v){ // row scalar multiplication
    vector<double> temp;
    for (auto &x : v)
        temp.push_back(k*x);
    return temp;
}
ostream& operator<<(ostream& os, const vector<double> &v){ // print a row
    if (not v.empty()) os << v[0];
    for (int i=1;i<v.size();++i){
        os << ' ' << v[i];
    }
    return os;
}
ostream& operator<<(ostream& os, const vector<vector<double>> &m){ // print a matrix
    for (auto &v: m)
        os << v << endl;
    return os;
}
struct ERO{
    int type;      // 1 for row swap, 2 for row addition
    int R1, R2;    // row index for the operation
    double k1, k2; // k1*R1 + k2*R2 -> R1
}; // if type 1 then k1, k2 are not used
ostream& operator<<(ostream& os, const ERO &ero){ // print ero steps
    if (ero.type == 1)
        os << "Swap rows " << ero.R1+1 << " and " << ero.R2+1;
    else if (ero.k2 == 0)
        os << "Multiply row " << ero.R1+1 << " by " << ero.k1;
    else os << "Add " << ero.k2 << " times row " << ero.R2+1 << " to row " << ero.R1+1;
    return os;
}
vector<vector<double>> transpose(const vector<vector<double>> &m){
    int M = m.size(), N = m.empty() ? 0 : m[0].size();
    vector temp(N, vector<double>(M, 0));
    for (int i=0;i<M;++i)
        for (int j=0;j<N;++j)
            temp[j][i] = m[i][j];
    return temp;
}
bool double_equal(const double &a, const double &b){ 
    return abs(a - b) < 1e-9;
}
bool operator==(const vector<double> &a, const vector<double> &b){
    if (a.size() != b.size()) return false;
    for (int i=0;i<a.size();++i)
        if (not double_equal(a[i], b[i])) return false;
    return true;
}
class matrix{
protected:
    int M, N, rk; // M rows, N columns
    vector<vector<double>> org; // original matrix
    vector<vector<double>> rref;
    vector<int> is_pivot; // is_pivot[i] = 1 if row i is a pivot column, 0 otherwise
    vector<vector<double>> basic_sol; // basis for the null space
    vector<ERO> gaussian; // gaussian elimination operations
public:
    int rows(){
        return M;
    }
    int cols(){
        return N;
    }
    int rank(){
        return rk;
    }
    int nullity(){
        return N - rk;
    }
    vector<vector<double>> original(){
        return org;
    }
    vector<vector<double>> RREF(){
        return rref;
    }
    void rowspace(){
        cout << "Dimension of row space: " << rk << endl;
        if (rk == 0){
            cout << "The row space has only the row vector 0" << endl;
            return;
        }
        cout << "A row space basis:" << endl;
        for (int i=0;i<rk;++i)
            cout << '[' << rref[i] << ']' << endl;
        cout << endl;
    }
    void colspace(){
        cout << "Dimension of column space: " << rk << endl;
        if (rk == 0){
            cout << "The column space has only the vector 0" << endl;
            return;
        }
        else if (rk == M)
            cout << "The column space is the whole space R^" << M << endl;
        cout << "A column space basis:" << endl;
        auto temp = transpose(org);
        for (int i=0;i<N;++i)
            if (is_pivot[i] != 0)
                cout << '[' << temp[i] << "]^T" << endl;
        cout << endl;
    }
    void nullspace(){
        cout << "Dimension of null space: " << N - rk << endl;
        if (rk == N){
            cout << "The null space has only the vector 0" << endl;
            return;
        }
        else if (rk == 0)
            cout << "The column space is the whole space R^" << N << endl;
        cout << "A null space basis:" << endl;
        for (auto &x : basic_sol)
            cout << '[' << x << "]^T" << endl;
        cout << endl;
    }
    void eval_gaussian(){ // follow wikipedia
        gaussian.clear();
        rref = org;
        int r = 0, c = 0; // pivot row and column
        while (r < M and c < N){
            int i = r;
            while (i < M and double_equal(rref[i][c], 0)) i++; // find pivot
            if (i == M) { // whole column is 0
                c++;
                continue;
            }
            if (i != r){
                gaussian.push_back(ERO{1, i, r, 0, 0});
                swap(rref[i], rref[r]);
                //print_rref();
            }
            if (not double_equal(rref[r][c], 1)){ // rescale ith row to leading 1
                gaussian.push_back(ERO{2, r, r, 1/rref[r][c], 0});
                rref[r] = (1/rref[r][c]) * rref[r];
                //print_rref();
            }
            for (i = 0; i < M; ++i){
                if (i == r or double_equal(rref[i][c], 0)) continue;
                gaussian.push_back(ERO{2, i, r, 1, -rref[i][c]}); 
                rref[i] = rref[i] + (-rref[i][c])*rref[r]; // ith row -= a[i][c] * pivot row
                //print_rref();
            }
            r++, c++;
        }
    }
    void eval_rank(){
        vector<double> zero(N, 0);
        for (int i=0;i<min(M, N);i++){
            if (rref[i] == zero){
                rk = i;
                return;
            }
        }
        rk = min(M, N);
    }
    void eval_pivot(){
        is_pivot.resize(M, 0);
        basic_sol.clear();
        vector<int> basic_var, free_var;
        for (int i=0;i<M;++i){
            for (int j=0;j<N;++j){
                if (not double_equal(rref[i][j], 0)){
                    is_pivot[j] = 1; // column j is a pivot column
                    break;
                }
            }
        }
        for (int i=0;i<N;++i){
            if (is_pivot[i] == 1) basic_var.push_back(i);
            else free_var.push_back(i);
        }
        for (auto &x : free_var){
            vector<double> sol(N, 0);
            sol[x] = 1;
            for (int i=0;i<basic_var.size();++i)
                sol[basic_var[i]] = -rref[i][x];
            basic_sol.push_back(sol);
        }
    }
    void solve(vector<double> b){ // solve Ax = b
        for (auto &ero : gaussian){
            if (ero.type == 1)
                swap(b[ero.R1], b[ero.R2]);
            else b[ero.R1] = ero.k1*b[ero.R1] + ero.k2*b[ero.R2];
        }
        vector<double> zero(N, 0);
        for (int i = 0; i < M; ++i) {
            if (rref[i] == zero and (not double_equal(b[i], 0))) {
                cout << "No solution" << endl;
                return;
            }
        }
        if (rk == N){
            cout << "Unique solution" << endl;
            vector<double> x;
            for (int i=0;i<N;++i)
                x.push_back(b[i]);
            cout << "x = [" << x << "]^T" << endl;
        }
        else{
            cout << "Infinite solutions" << endl;
            cout << "General solution: x = x_p";
            for (int i=0;i<basic_sol.size();++i)
                cout << " + t_" << i+1 << "*x_" << i+1;
            cout << endl; 
            cout << "Where the t\'s are real numbers and" << endl;
            vector<double> x_p(N, 0);
            for (int i=0, j=0;i<N;++i){
                if (is_pivot[i] == 1)
                    x_p[i] = b[j++];
            }

            cout << "x_p = [" << x_p << "]^T" << endl;
            for (int i=0;i<basic_sol.size();++i){
                cout << "x_" << i+1 << " = [" << basic_sol[i] << "]^T" << endl;
            }
            cout << endl;
        }
    }
    matrix(const vector<vector<double>>& mat)
    : M(mat.size()), N(mat.empty() ? 0 : mat[0].size()), org(mat) {
        eval_gaussian();
        eval_rank();
        eval_pivot();
    }
};
class sqmatrix :  public matrix{
private:
    double det;
    vector<vector<double>> inv; // empty if not invertible
public:
    double determinant(){
        return det;
    }
    vector<vector<double>> inverse(){
        return inv;
    }
    void eval_det(){
        if (rk != N){
            det = 0;
            return;
        }
        det = 1;
        for (auto &ero : gaussian){
            if (ero.type == 1){ // row swap
                det *= -1;
            }
            else{ // k1*R1 + k2*R2 -> R1
                det *= (1/ero.k1);
            }
        }
    }
    void eval_inv(){
        if (rk != N) {
            inv.clear();
            return;
        }
        inv = identity_matrix(N);
        for (auto &ero : gaussian){
            if (ero.type == 1)
                swap(inv[ero.R1], inv[ero.R2]);
            else inv[ero.R1] = ero.k1*inv[ero.R1] + ero.k2*inv[ero.R2];
        }
    }
    sqmatrix (const vector<vector<double>> &mat): matrix(mat){
        if (M != N) return;
        eval_det();
        eval_inv();
    }
};
int get_pos_int(const string &prompt){ // cin >> positive integer
    string input, token;
    while (true){
        cout << prompt;
        cout << "(press c to cancel) ";
        getline(cin, input);
        if (input == "c") return -1;
        if (input.empty()) continue;
        bool flag = true;
        for (char c : input){
            if (c < '0' or c > '9') {
                cout << "E: Invalid input" << endl;
                flag = false;
                break;
            }
        }
        if (flag == false) continue;
        try {
            int n = stoi(input);
            if (n <= 0) cout << "E: Input is not positive" << endl;
            else if (n > 1000) cout << "E: Input out of range" << endl;
            else return n;
        }
        catch (const invalid_argument &e) {
            cout << "E: Invalid input" << endl;
        }
        catch (const out_of_range &e) {
            cout << "E: Input out of range" << endl;
        }
    }
}
vector<double> get_row(int M, const string &prompt){ // cin >> row of M entries
    string input, token;
    vector<double> temp;
    while (true){
        bool flag = true;
        temp.clear();
        cout << prompt;
        cout << "(press c to cancel) ";
        getline(cin, input);
        if (input == "c") return temp; // cancellation
        stringstream ss(input);
        while (ss >> token) {
            try {
                double x = stod(token);
                if (not isfinite(x)){
                    cout << "E: Invalid input" << endl;
                    flag = false;
                    break;
                }
                temp.push_back(stod(token));
            }
            catch (const invalid_argument &e) {
                cout << "E: Invalid input" << endl;
                flag = false;
                break;
            }
            catch (const out_of_range &e) {
                cout << "E: Input out of range" << endl;
                flag = false;
                break;
            }
        }
        if (flag == false) continue;
        else if (temp.size() != M) {
            cout << "E: Number of entries should be " << M << endl;
            continue;
        }
        else return temp;
    }
}
matrix get_matrix(){
    vector<vector<double>> mat;
    int rows = get_pos_int("Enter the number of rows: ");
    if (rows == -1) return matrix(mat); // cancellation
    int cols = get_pos_int("Enter the number of columns: ");
    if (cols == -1) return matrix(mat); // cancellation
    for (int i = 0; i < rows; ++i){
        auto temp = get_row(cols, "Enter row " + to_string(i+1) + ": ");
        if (temp.empty()){ // cancellation
            mat.clear();
            return matrix(mat);
        }
        else mat.push_back(temp);
    }
    return matrix(mat);
}
void print_help(){
    cout << "Available operations:" << endl;
    cout << "h: print this help message" << endl;
    cout << "q: quit" << endl;
    cout << "d: delete this matrix" << endl;
    cout << "org: print the original matrix" << endl;
    cout << "rref: print the reduced row echelon form of the matrix" << endl;
    cout << "inv: print the inverse matrix" << endl;
    cout << "size: matrix size" << endl;
    cout << "rank: matrix rank" << endl;
    cout << "nullity: matrix nullity" << endl;
    cout << "det: matrix determinant" << endl;
    cout << "rowspace: matrix row space" << endl;
    cout << "colspace: matrix column space" << endl;
    cout << "nullspace: matrix null space" << endl;
    cout << "solve: solve the linear system Ax=b" << endl;
    cout << endl;
}
int main(){
    cout << fixed << setprecision(4);
    cout << "Welcome to the Matrix Calculator" << endl;
    string input;
    while (true){
        cout << "Press n to enter a new matrix, q to exit: " << endl;
        getline(cin, input);
        if (input == "q"){ // quit
            cout << "Exit" << endl;
            return 0;
        } 
        else if (input != "n"){ // invalid input
            cout << "E: Invalid input" << endl;
            continue;
        }
        else{
            matrix mat = get_matrix();
            if (mat.original().empty()) continue;
            bool is_square = (mat.rows() == mat.cols());
            sqmatrix sqmat(mat.original());
            cout << endl << "Matrix created" << endl;
            print_help();
            while (true) {
                cout << "Enter an operation: (press h for help) ";
                getline(cin, input);
                if (input == "q"){ // quit
                    cout << "Exit" << endl;
                    return 0;
                }
                else if (input == "h"){ // print help
                    print_help();
                }
                else if (input == "d"){ // delete matrix
                    cout << "Matrix deleted." << endl;
                    break;
                } 
                else if (input == "org") { // print original matrix
                    cout << mat.original() << endl;
                }
                else if (input == "rref") { // print RREF
                    cout << mat.RREF() << endl;
                    
                }
                else if (input == "inv") { // print inverse matrix
                    if (not is_square) {
                        cout << "E: Matrix is not square." << endl;
                        continue;
                    }
                    if (sqmat.determinant() == 0)
                        cout << "E: Matrix is not invertible." << endl;
                    else cout << sqmat.inverse() << endl;
                }
                else if (input == "size") { // print size
                    cout << "Size: " << mat.rows() << 'x' << mat.cols() << endl;
                } 
                else if (input == "rank") { // print rank
                    cout << "Rank: " << mat.rank() << endl;
                }
                else if (input == "nullity") { // print nullity
                    cout << "Nullity: " << mat.nullity() << endl;
                }
                else if (input == "det") { // print determinant
                    if (not is_square) {
                        cout << "E: Matrix is not square." << endl;
                        continue;
                    }
                    else cout << "Determinant: " << sqmat.determinant() << endl;
                }
                else if (input == "rowspace") { // print row space
                    mat.rowspace();
                } 
                else if (input == "colspace") { // print column space
                    mat.colspace();
                } 
                else if (input == "nullspace") { // print null space
                    mat.nullspace();
                }
                else if (input == "solve") { // solve linear system Ax=b
                    cout << "Solve the linear system Ax=b" << endl;
                    vector<double> b = get_row(mat.rows(), "Enter the right-hand side vector b: ");
                    if (b.empty()) continue;
                    mat.solve(b);
                }
                else{ // invalid input
                    cout << "E: Invalid input" << endl;
                }
            }
        }
    }
    return 0;
}