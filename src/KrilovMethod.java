import java.util.ArrayList;
import java.util.List;

public class KrilovMethod {
    private double vectorQ[];
    private double matrC[][];
    private double vectorC[];
    private double matrA[][];
    private List<double[]> eigenVectors;
    private List<double[]> discrepancy;
    private double lambda[] = {0.191009010539547, 0.3835583164926, 0.597703453394, 0.879661472553, 1.1447561970210};
    private int n;

    KrilovMethod(double matrA[][]){
        n = matrA.length;
        this.matrA = new double[n][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                this.matrA[i][j] = 0;
                for (int k = 0; k < n; k++)
                    this.matrA[i][j] += matrA[k][i] * matrA[k][j];
            }
        }
    }


    private double [] multiply(double[][] matr, double[] vector){
        double result[] = new double[n];
        for(int i = 0; i < n; i++)
        {
            result[i]=0;
            for(int j = 0; j < n; j++)
            {
                result[i] += matr[i][j]*vector[j];
            }
        }
        return result;
    }
    private double[][] multiply(double[][] matrA, double[][] matrB){
        double result[][] = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = 0;
                for (int k = 0; k < n; k++)
                    result[i][j] += matrA[i][k] * matrB[k][j];
                }
            }
        return result;
    }

    private double[] multiply(double[] vector, double lambda){
        double[] temp = new double[n];
        for(int i = 0; i < n; i++){
            temp[i] = vector[i] * lambda;
        }
        return temp;
    }
    private double[] minus(double[] vectorA, double[] vectorB){
        double[] temp = new double[n];
        for(int i = 0; i < n; i++){
            temp[i] = vectorA[i] - vectorB[i];
        }
        return temp;
    }
    private double[][] pow(double[][] matr, int pow){
        if(pow == 1){
            return matr;
        }
        double resMatr[][] = new double[n][n];
        resMatr = multiply(matr, matr);
        for(int count = 2; count < pow; count++){
            resMatr = multiply(resMatr, matr);
        }
        return resMatr;
    }

    private void createMatrC(){
        matrC = new double[n][n];
        double vectC0[] = new double[n];
        vectC0[0] = 1;
        for(int i = 1; i < n; i++){
            vectC0[i] = 0;
        }
        for(int j = 0; j < n; j++){
            for(int i = 0; i < n; i++){
                if(j == n-1){
                    matrC[i][j] = vectC0[i];
                }
                else {
                    matrC[i][j] = multiply(pow(matrA, n - 1 - j), vectC0)[i];
                }
            }
        }
        vectorC = new double[n];
        for(int i = 0; i < n; i++){
            vectorC[i] = multiply(pow(matrA, n), vectC0)[i];
        }
    }


    public void krilovMethod(){
        createMatrC();
        GaussMethod gs = new GaussMethod(matrC.clone(), vectorC);
        gs.strightStep();
        gs.createVectorX();
        vectorQ = new double[n];
        System.out.println("Q:");
        for(int i = 0; i < n; i++) {
            vectorQ[i] = gs.getVectorX()[i];
            System.out.println(vectorQ[i]);
        }
        System.out.println();
    }

    private double[] createEigenVector(double lambda){
        double beta[] = new double[n];
        beta[0] = 1;
        for(int j = 1; j < n; j++){
            beta[j] = lambda * beta[j - 1] - vectorQ[j - 1];
        }
        double[][] CT = transp(matrC);
        double[] eigen = new double[n];
        for (int i = 0; i < n; ++i) {
            eigen = plus(eigen, multiply(CT[i], beta[i]));
        }
        return eigen;
    }

    public void createEigenVectors(){
        eigenVectors = new ArrayList<>();
        for(int i = 0; i < n; i++){
            eigenVectors.add(createEigenVector(lambda[i]));
        }
    }
    public void printEigenVectors(){
        int count = 0;
        for(double[] vector: eigenVectors){
            System.out.println("Eigen vector for lambda = " + lambda[count]);
            for(double item : vector){
                System.out.format("%5f  ", item);
            }
            System.out.println();
            count++;
        }
    }

    public void printDiscrepancy(){
        double [] polinom = new double[n + 1];
        polinom[0] = 1;
        for(int i = 1; i < n+1; i++){
            polinom[i] = -vectorQ[i-1];
        }
        double discrepancy[] = new double[n];
        for(int i = 0; i < n; i++){
            discrepancy[i] = 0;
            for(int j = 0; j < n; j++){
                discrepancy[i] += Math.pow(lambda[i], n-j) * polinom[j];
            }
            discrepancy[i] += polinom[n];
            System.out.println("discrepancy for lambda = " + lambda[i]);
            System.out.format("%25s", discrepancy[i] + "\n");
        }
    }

    private double[] createVectorDiscrepancy(double[] eigenVector, double lambda){
        return minus(multiply(matrA, eigenVector), multiply(eigenVector, lambda));
    }

    public void createVectorsDiscrepancy(){
        discrepancy = new ArrayList<>();
        for(int i = 0; i < n; i++){
            discrepancy.add(createVectorDiscrepancy(eigenVectors.get(i), lambda[i]));
        }
    }

    public void printVectorsDiscrepancy(){
        int count = 0;
        System.out.println("\nEigen vectors discrepancy\n");
        for(double[] vector: discrepancy){
            System.out.println("discrepancy for lambda = " + lambda[count]);
            for(double item : vector){
                System.out.format("%e ", item);
            }
            System.out.println();
            count++;
        }
    }

    public double[][] transp(double[][] matrix) {
        double[][] result = new double[n][n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result[i][j] = matrix[j][i];
            }
        }
        return result;
    }

    public double[] plus(double[] v1, double[] v2) {
        double[] result = new double[n];
        if (v1.length == v2.length) {
            for (int i = 0; i < n; ++i) {
                result[i] = v1[i] + v2[i];
            }
        }
        return result;
    }
}
