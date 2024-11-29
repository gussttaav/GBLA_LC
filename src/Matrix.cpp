#include "Matrix.h"

namespace GBLA_LC
{

    Matrix::Matrix()
    {
        cols = 0;
        rows = 0;

        m = NULL;
    }

    Matrix::Matrix(int r, int c)
    {
        rows = r;
        cols = c;

        m = new GF*[rows];
        for(int i = 0; i < rows; i++)
            m[i] = new GF[cols];
    }

    Matrix::Matrix(int r, int c, GF &e)
    {
        rows = r;
        cols = c;

        m = new GF*[rows];
        for(int i = 0; i < rows; i++)
            m[i] = new GF[cols];

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                m[i][j] = e;
    }

    Matrix::Matrix(const Matrix &mat)
    {
        rows = mat.rows;
        cols = mat.cols;

        m = new GF*[rows];
        for(int i = 0; i < rows; i++)
            m[i] = new GF[cols];

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                m[i][j] = mat.m[i][j];
    }

    Matrix& Matrix::operator=(const Matrix& mat)
    {
        if(this == &mat)
            return *this;

        if(m)
        {
            for(int i = 0; i < rows; i++)
                delete [] m[i];
            delete [] m;
        }

        m = new GF*[mat.rows];
        for(int i = 0; i < mat.rows; i++)
            m[i] = new GF[mat.cols];

        rows = mat.rows;
        cols = mat.cols;

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                m[i][j] = mat.m[i][j];

        return *this;
    }

    Matrix::~Matrix()
    {
        if(m)
        {
            for(int i = 0; i < rows; i++)
                delete [] m[i];
            delete [] m;
        }
    }


    void Matrix::swap_rows(int row1, int row2)
    {
        GF tmp;
        for(int i = 0; i < cols; i++)
        {
            tmp = m[row1][i];
            m[row1][i] = m[row2][i];
            m[row2][i] = tmp;
        }
    }

    void Matrix::swap_cols(int col1, int col2)
    {
        GF tmp;
        for(int i = 0; i < rows; i++)
        {
            tmp = m[i][col1];
            m[i][col1] = m[i][col2];
            m[i][col2] = tmp;
        }
    }


    void Matrix::permute_cols(const std::vector<int> &perm)
    {
        Matrix t = *this;
        for(int j = 0; j < cols; j++)
            if(j != perm[j])
            {
                for(int i = 0; i < rows; i++)
                    m[i][perm[j]] = t(i, j);
            }
    }

    int Matrix::rank() const
    {
        GF mult;
        int rnk = 0;
        Matrix B(*this);
        int rp = 0, cp = 0;

        for(; (rp < rows-1) && (cp < cols); rp++, cp++)
        {
            if(B(rp,cp) != 0)
            {
                rnk++;
                for(int k = rp+1; k < rows; k++)
                {
                    if(B(k,cp) != 0)
                    {
                        mult = B(k,cp) / B(rp,cp);
                        for(int j = cp; j < cols; j++)
                            B(k,j) -= mult * B(rp,j);
                    }
                }
            }
            else{
                bool swapped = false;

                for(int k = rp+1; k < rows; k++)
                {
                    if(B(k,cp) != 0)
                    {
                        swapped = true;
                        B.swap_rows(rp, k);
                        break;
                    }
                }

                if(swapped)
                    cp--;
                rp--;
            }
        }

        //Analyze last row
        for(int j = cols-1; j >= 0; j--)
            if(B(rp,j) != 0)
               return ++rnk;

        return rnk;
    }

    std::vector<GF> operator*(const std::vector<GF> &v, const Matrix &mat)
    {
        assert(v.size() == (unsigned)mat.rows);

        std::vector<GF> result(mat.cols, 0);

        for(int i = 0; i < mat.cols; i++)
            for(int j = 0; j < mat.rows; j++)
                result[i] += v[j]*mat.m[j][i];

        return result;
    }

    std::ostream& Matrix::output(std::ostream& out) const
    {
        int max_width = 0, i, j, k, l, t = cols-1;
        std::vector<std::string> vec(rows*cols);

        for(i = 0, k = 0; i < rows; i++)
        {
            for(j = 0; j < cols; j++, k++)
            {
                std::ostringstream os;
                os << m[i][j];
                vec[k] = os.str();

                if(max_width < (int)vec[k].length())
                    max_width = vec[k].length();
            }
        }

        for(i = 0, k = 0; i < rows; i++)
        {
            out << "[";
            for(j = 0; j < t; j++, k++)
            {
                // add spaces around the string to center it
                l = max_width - vec[k].length();
                if(l % 2)
                    vec[k] = " " + vec[k];

                for(l = l/2; l > 0; l--)
                    vec[k] = " " + vec[k] + " ";
                // output the centered string
                out << vec[k] << " ";
            }
            // add spaces around the string to center it
            l = max_width - vec[k].length();
            if(l % 2)
                vec[k] = " " + vec[k];

            for(l = l/2; l > 0; l--)
                vec[k] = " " + vec[k] + " ";

            // output the centered string
            out << vec[k++] << "]" << std::endl;
        }

        return out;
    }

    std::ostream &operator<<(std::ostream &out, const Matrix &mat)
    {   return mat.output(out); }

}
