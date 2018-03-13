import java.io.IOException;
import java.util.function.DoubleUnaryOperator;

/**
 * Created by Миша on 04.03.2018.
 */
public class Main {
    public static void main(String[] args) throws Exception {
        int n = 3;
        int N = 13; //accuracy of taking integrals
        double a = -0.7978525, b = -0.5973559;
//        double a = -1, b = 1;
//        double k = 1.5919,        l = 5.4236;
        double k = 15.7503, l = 19.58201;
        DoubleUnaryOperator p = x -> k * x + l;
        DoubleUnaryOperator q = x -> k * k * (1 / (k * x + l) - k * x);
        double p_min = p.applyAsDouble(a);
        double p_max = p.applyAsDouble(b);
        double q_min = q.applyAsDouble(b);
        double q_max = q.applyAsDouble(a);
        double nu_1 = Math.PI / 2;
        double nu_2 = Math.PI;
        double lambda1_min = 4 * nu_1 * nu_1 / Math.pow(b - a, 2) * p_min + q_min;
        double lambda1_max = 4 * nu_1 * nu_1 / Math.pow(b - a, 2) * p_max + q_max;
        double lambda2_min = 4 * nu_2 * nu_2 / Math.pow(b - a, 2) * p_min + q_min;
        double lambda2_max = 4 * nu_2 * nu_2 / Math.pow(b - a, 2) * p_max + q_max;
        System.out.println(lambda1_min);
        System.out.println(lambda1_max);
//        System.out.println(lambda2_min);
//        System.out.println(lambda2_max);
        double c_1 = Math.sqrt(2 / (b - a));
        double c_2 = c_1;
        DoubleUnaryOperator y_1 = t -> c_1 * Math.cos(2 * nu_1 * t / (b - a) - (a + b) * nu_1 / (b - a));
        DoubleUnaryOperator y_2 = t -> c_2 * Math.sin(2 * nu_2 * t / (b - a) - (a + b) * nu_2 / (b - a));

        DoubleUnaryOperator L_1_min = t -> -(p_min * (-c_1 * 4 * nu_1 * nu_1 / Math.pow(b - a, 2) * Math.cos(2 * nu_1 * t / (b - a) - (a + b) * nu_1 / (b - a)))) + q_min * c_1 * Math.cos(2 * nu_1 * t / (b - a) - (a + b) * nu_1 / (b - a));
        DoubleUnaryOperator L_1_max = t -> -(p_max * (-c_1 * 4 * nu_1 * nu_1 / Math.pow(b - a, 2) * Math.cos(2 * nu_1 * t / (b - a) - (a + b) * nu_1 / (b - a)))) + q_max * c_1 * Math.cos(2 * nu_1 * t / (b - a) - (a + b) * nu_1 / (b - a));

        DoubleUnaryOperator L_2_min = t -> -(p_min * (-c_2 * 4 * nu_2 * nu_2 / Math.pow(b - a, 2) * Math.sin(2 * nu_2 * t / (b - a) - (a + b) * nu_2 / (b - a)))) + q_min * c_2 * Math.sin(2 * nu_2 * t / (b - a) - (a + b) * nu_2 / (b - a));
        DoubleUnaryOperator L_2_max = t -> -(p_max * (-c_2 * 4 * nu_2 * nu_2 / Math.pow(b - a, 2) * Math.sin(2 * nu_2 * t / (b - a) - (a + b) * nu_2 / (b - a)))) + q_max * c_2 * Math.sin(2 * nu_2 * t / (b - a) - (a + b) * nu_2 / (b - a));

//        System.out.println(L_1_min.applyAsDouble((a+b)/2)-lambda1_min*y_1.applyAsDouble((a+b)/2));
//        System.out.println(L_1_max.applyAsDouble((a+b)/2)-lambda1_max*y_1.applyAsDouble((a+b)/2));
//        System.out.println(L_2_min.applyAsDouble((a+b)/2)-lambda2_min*y_2.applyAsDouble((a+b)/2));
//        System.out.println(L_2_max.applyAsDouble((a+b)/2)-lambda2_max*y_2.applyAsDouble((a+b)/2));

//        System.out.println(y_1.applyAsDouble(a));
//        System.out.println(y_1.applyAsDouble(b));
//        System.out.println(y_2.applyAsDouble(a));
//        System.out.println(y_2.applyAsDouble(b));


        //part 2:
        Gauss g = new Gauss();
        DoubleUnaryOperator energ1 = x -> p.applyAsDouble(x) * Math.pow(-c_1 * 2 * nu_1 / (b - a) * Math.cos(2 * nu_1 * x / (b - a) - (a + b) * nu_1 / (b - a)), 2) + q.applyAsDouble(x) * Math.pow(c_1 * Math.cos(2 * nu_1 * x / (b - a) - (a + b) * nu_1 / (b - a)), 2);
        DoubleUnaryOperator scalar1 = x -> Math.pow(c_1 * Math.cos(2 * nu_1 * x / (b - a) - (a + b) * nu_1 / (b - a)), 2);
        double numerator1 = g.takeIntegral(a, b, energ1, N);
        double denominator1 = g.takeIntegral(a, b, scalar1, N);
        double lambda_1 = numerator1 / denominator1;

//        System.out.println(lambda_1);

        DoubleUnaryOperator energ2 = x -> p.applyAsDouble(x) * Math.pow(c_2 * 2 * nu_2 / (b - a) * Math.sin(2 * nu_2 * x / (b - a) - (a + b) * nu_2 / (b - a)), 2) + q.applyAsDouble(x) * Math.pow(c_2 * Math.sin(2 * nu_2 * x / (b - a) - (a + b) * nu_2 / (b - a)), 2);
        DoubleUnaryOperator scalar2 = x -> Math.pow(c_2 * Math.sin(2 * nu_2 * x / (b - a) - (a + b) * nu_2 / (b - a)), 2);
        double numerator2 = g.takeIntegral(a, b, energ2, N);
        double denominator2 = g.takeIntegral(a, b, scalar2, N);
        double lambda_2 = numerator2 / denominator2;

//part 3
        DoubleUnaryOperator[] P_2_2 = new DoubleUnaryOperator[n];
        P_2_2[0] = x -> 1;
        P_2_2[1] = x -> 3 * x;
        for (int i = 0; i < n - 2; i++) {
            double j = i;
            int finalI = i;
            P_2_2[i + 2] = x -> (j + 4) * ((2 * j + 7) * x * P_2_2[finalI + 1].applyAsDouble(x) - (j + 3) * P_2_2[finalI].applyAsDouble(x)) / (j + 2) / (j + 6);
        }
        DoubleUnaryOperator[] P_2_2_norm = new DoubleUnaryOperator[n];
        for (int i = 0; i < n; i++) {
            double j = i;
            int finalI = i;
            P_2_2_norm[i] = t -> d(j) * P_2_2[finalI].applyAsDouble(t);
        }

        DoubleUnaryOperator[] omega = new DoubleUnaryOperator[n];
        for (int i = 0; i < n; i++) {
            int finalI = i;
            omega[i] = x -> (1 - x * x) * P_2_2_norm[finalI].applyAsDouble(x);
        }

        DoubleUnaryOperator[] P_3_3 = new DoubleUnaryOperator[n];
        P_3_3[0] = x -> 1;
        P_3_3[1] = x -> 4 * x;
        for (int i = 0; i < n - 2; i++) {
            double j = i;
            int finalI = i;
            P_3_3[i + 2] = x -> (j + 5) * ((2 * j + 9) * x * P_3_3[finalI + 1].applyAsDouble(x) - (j + 4) * P_3_3[finalI].applyAsDouble(x)) / (j + 2) / (j + 8);
        }
        DoubleUnaryOperator[] dP_2_2 = new DoubleUnaryOperator[n];
        dP_2_2[0] = x -> 0;
        for (int i = 1; i < n; i++) {
            double j = i;
            int finalI = i;
            dP_2_2[i] = x -> (j + 5) / 2 * P_3_3[finalI - 1].applyAsDouble(x);
        }

        DoubleUnaryOperator[] dotOmega = new DoubleUnaryOperator[n];
        for (int i = 0; i < n; i++) {
            double j = i;
            int finalI = i;
            dotOmega[i] = x -> (-2 * x) * P_2_2_norm[finalI].applyAsDouble(x) + (1 - x * x) * d(j) * dP_2_2[finalI].applyAsDouble(x);
        }


        double[][] G_L = new double[n][n];
        double[][] G = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i + 1; j++) {
                int finalI = i;
                int finalJ = j;
                DoubleUnaryOperator energ = t -> (p.applyAsDouble(t) * dotOmega[finalI].applyAsDouble(2 * t / (b - a) - (a + b) / (b - a)) * dotOmega[finalJ].applyAsDouble(2 * t / (b - a) - (a + b) / (b - a)) + q.applyAsDouble(t) * omega[finalI].applyAsDouble(2 * t / (b - a) - (a + b) / (b - a)) * omega[finalJ].applyAsDouble(2 * t / (b - a) - (a + b) / (b - a)));
                DoubleUnaryOperator scalar = t -> 2 / (b - a) * omega[finalI].applyAsDouble(2 * t / (b - a) - (a + b) / (b - a)) * omega[finalJ].applyAsDouble(2 * t / (b - a) - (a + b) / (b - a));
                G_L[i][j] = g.takeIntegral(a, b, energ, N);
//                G[i][j] = g.takeIntegral(a, b, scalar, N);
                if (i == j)
                    G[i][j] = 1;
                else G[i][j] = 0;
            }
        }
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                G_L[i][j] = G_L[j][i];
                G[i][j] = G[j][i];
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 5; j++) {
                System.out.print(dotOmega[i].applyAsDouble(a + j / 4.0 * (b - a)) + "   ");
            }
            System.out.println();
        }

        Matrix matrix = new Matrix(G_L);
//        matrix.printf();
        ScalarProductMethod method = new ScalarProductMethod();
        System.out.println(1 / method.getMaxEigenvalue(matrix.getInverse()));

        System.out.println("answer = " + k * k * l);
    }


    private static double d(double j) {
        return 0.25 * Math.sqrt((j + 3) * (j + 4) * (2 * j + 5) / (2 * (j + 1) * (j + 2)));
    }
}
