#ifndef MATRIX_H
#define MATRIX_H

#include "GF.h"

using namespace GaloisField;

namespace GBLA_LC
{

    class Matrix
    {
        public:
            Matrix();
            Matrix(int r, int c);
            Matrix(int r, int c, GF &e);
            Matrix(const Matrix &mat);
            ~Matrix();

            Matrix& operator=(const Matrix& mat);
            GF& operator()(int i, int j) { return m[i][j]; }
            const GF& operator()(int i, int j) const { return m[i][j]; }
            friend std::vector<GF> operator*(const std::vector<GF> &v, const Matrix &mat);

            int get_rows() const { return rows; }
            int get_cols() const { return cols; }
            void permute_cols(const std::vector<int> &perm);

            int rank() const;
            bool is_empty() const { return rows == 0 && cols == 0; }

            std::ostream &output(std::ostream& out) const;
            friend std::ostream &operator<<(std::ostream &out, const Matrix& mat);

        private:
            int rows, cols;
            GF **m;

            void swap_rows(int row1, int row2);
            void swap_cols(int col1, int col2);
    };

}

#endif // MATRIX_H
