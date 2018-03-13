
/**
 * Created by Миша on 26.09.2017.
 */

public class Kramer {


    public double[] solve(Matrix A, double[] b) throws Exception {
        int n = A.matrix.length;

        Matrix[] masOfMatrix = new Matrix[n];
        for (int i = 0; i < n; i++) {
            masOfMatrix[i] = new Matrix(n);
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    masOfMatrix[i].matrix[j][k] = A.matrix[j][k];
                }
            }
            for (int j = 0; j < n; j++) {
                masOfMatrix[i].matrix[j][i] = b[j];
            }
        }


        double[] x = new double[n];
        double d = A.det();
        for (int i = 0; i < n; i++) {
            x[i] = masOfMatrix[i].det() / d;
        }


        return x;
    }

}
