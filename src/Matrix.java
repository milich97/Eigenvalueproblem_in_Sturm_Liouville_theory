
/**
 * Created by Миша on 26.09.2017.
 */
public class Matrix {
    public double[][] matrix;

    public Matrix(int n) {
        matrix = new double[n][n];
    }

    public Matrix(double[][] A) {
        matrix = new double[A.length][A.length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A.length; j++) {
                matrix[i][j] = A[i][j];
            }
        }

    }

    public Matrix() {
    }

    public Matrix plus(Matrix m) {
        Matrix res = new Matrix(m.matrix.length);
        for (int i = 0; i < m.matrix.length; i++) {
            for (int j = 0; j < m.matrix.length; j++) {
                res.matrix[i][j] = this.matrix[i][j] + m.matrix[i][j];
            }
        }
        return res;
    }

    public Matrix getInverse() throws Exception {
        Matrix inverseMatrix = new Matrix(matrix.length);
        for (int i = 0; i < matrix.length; i++) {
            Vector b = new Vector(matrix.length);
            for (int j = 0; j < matrix.length; j++) {
                b.vector[j] = kroneker(i, j);
            }
            Kramer kramer = new Kramer();
            b.vector = kramer.solve(this, b.vector);
            for (int j = 0; j < b.vector.length; j++) {
                inverseMatrix.matrix[j][i] = b.vector[j];
            }
        }
        return inverseMatrix;
    }

    public Matrix mulToNumber(double a) {
        Matrix res = new Matrix(matrix.length);
        for (int i = 0; i < this.matrix.length; i++) {
            for (int j = 0; j < this.matrix.length; j++) {
                res.matrix[i][j] = a * this.matrix[i][j];
            }
        }
        return res;
    }

    private int kroneker(int i, int j) {
        if (i == j) return 1;
        return 0;
    }

    public Vector mulToVector(Vector x0) {
        Vector res = new Vector(x0.vector.length);
        for (int i = 0; i < x0.vector.length; i++) {
            for (int j = 0; j < x0.vector.length; j++) {
                res.vector[i] += this.matrix[i][j] * x0.vector[j];
            }

        }
        return res;
    }

    public Matrix mulToMatrix(Matrix A) {
        Matrix res = new Matrix(A.matrix.length);
        for (int i = 0; i < A.matrix.length; i++) {
            for (int j = 0; j < A.matrix.length; j++) {
                for (int k = 0; k < A.matrix.length; k++) {
                    res.matrix[i][j] += this.matrix[i][k] * A.matrix[k][j];
                }
            }

        }
        return res;
    }

    public double det() {
        double res;
        int N = matrix.length;
        if (N == 1)
            res = matrix[0][0];
        else if (N == 2)
            res = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
        else {
            res = 0;
            for (int j1 = 0; j1 < N; j1++) {
                Matrix m = new Matrix(N);
                m.matrix = generateSubArray(matrix, N, j1);
                res += Math.pow(-1.0, 1.0 + j1 + 1.0) * matrix[0][j1] * m.det();
            }
        }
        return res;
    }

    public double[][] generateSubArray(double A[][], int N, int j1) {
        double[][] m = new double[N - 1][];
        for (int k = 0; k < (N - 1); k++)
            m[k] = new double[N - 1];

        for (int i = 1; i < N; i++) {
            int j2 = 0;
            for (int j = 0; j < N; j++) {
                if (j == j1)
                    continue;
                m[i - 1][j2] = A[i][j];
                j2++;
            }
        }
        return m;
    }

    public Matrix minus(Matrix m) {
        Matrix res = new Matrix(m.matrix.length);
        for (int i = 0; i < m.matrix.length; i++) {
            for (int j = 0; j < m.matrix.length; j++) {
                res.matrix[i][j] = this.matrix[i][j] - m.matrix[i][j];
            }
        }
        return res;
    }

    public Matrix transpose() {
        Matrix m = new Matrix(matrix.length);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                m.matrix[i][j] = matrix[j][i];
            }
        }
        return m;
    }

    public void fillE() {
        for (int i = 0; i < this.matrix.length; i++) {
            for (int j = 0; j < this.matrix.length; j++) {
                if (i == j) {
                    this.matrix[i][j] = 1;
                } else {
                    this.matrix[i][j] = 0;
                }
            }
        }
    }

    public void printf() {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }
}
