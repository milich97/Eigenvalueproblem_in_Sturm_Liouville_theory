
import java.util.Random;

/**
 * Created by Миша on 26.09.2017.
 */
public class Vector {
    double[] vector;

    public Vector(int n) {
        vector = new double[n];
    }

    public Vector(double v1, double v2, double v3) {
        vector = new double[3];
        vector[0] = v1;
        vector[1] = v2;
        vector[2] = v3;
    }


    public void printf() {
        for (int i = 0; i < this.vector.length; i++) {
            System.out.println(this.vector[i]);
        }
    }

    public void normalize() {
        int i = 0;
        for (int j = 1; j < vector.length; j++) {
            if (Math.abs(vector[j]) > Math.abs(vector[i])) i = j;
        }
        double q = vector[i];
        for (int j = 0; j < vector.length; j++) {
            vector[j] /= q;
        }
    }

    public void mulToNumber(double alpha) {
        for (int i = 0; i < this.vector.length; i++) {
            this.vector[i] *= alpha;
        }
    }

    public Vector plus(Vector x0) {
        Vector res=new Vector(x0.vector.length);
        for (int i = 0; i <x0.vector.length ; i++) {
            res.vector[i]=this.vector[i]+x0.vector[i];
        }
        return res;
    }

    public Vector minus(Vector x0) {
        Vector res=new Vector(x0.vector.length);
        for (int i = 0; i <x0.vector.length ; i++) {
            res.vector[i]=this.vector[i]-x0.vector[i];
        }
        return res;
    }

    public double norm() {
        double r = 0;
        for (int i = 0; i < this.vector.length; i++) {
            r += Math.abs(this.vector[i]);
        }
        return r;
    }
}
