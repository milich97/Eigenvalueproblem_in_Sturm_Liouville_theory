
/**
 * Created by Миша on 04.11.2017.
 */
public class ScalarProductMethod {
    public double getMaxEigenvalue(Matrix A) {
        double eps = 0.001;
        Vector v0 = new Vector(A.matrix.length);
        for (int i = 0; i < v0.vector.length; i++) {
            v0.vector[i]=1;
        }
        Vector v1 = A.mulToVector(v0);
        double maxEigenvalue0 = Double.MAX_VALUE;
        double s1 = 0;
        double s2 = 0;
        for (int i = 0; i < v0.vector.length; i++) {
            s1 += v0.vector[i] * v1.vector[i];
            s2 += v0.vector[i] * v0.vector[i];
        }
        double maxEigenvalue1 = s1 / s2;
        int k=0;
        while (Math.abs(maxEigenvalue0 - maxEigenvalue1) > eps) {
            k++;
            s1 = 0;
            s2 = 0;
            maxEigenvalue0 = maxEigenvalue1;
            v0 = v1;
            v0.normalize();
            v1 = A.mulToVector(v0);
            for (int i = 0; i < v0.vector.length; i++) {
                s1 += v0.vector[i] * v1.vector[i];
                s2 += v0.vector[i] * v0.vector[i];
            }
            maxEigenvalue1 = s1 / s2;
        }
//        System.out.println("scalar.method "+k);
//        System.out.println(maxEigenvalue1);
        return maxEigenvalue1;
    }
}
