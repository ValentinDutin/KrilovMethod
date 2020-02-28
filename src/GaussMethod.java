import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.abs;

public class GaussMethod {
    private double[][] resultMatrA;
    private double[][] matrA;
    private double detA;
    private List<Double> mainElements;
    private int numberOfPermutations = 0;
    private double[] vectorX;
    private double[] vectorB;
    private double[] copyVectorB;
    private int n;
    private int[][] queue;
    private double[] discrepancy;
    private double[][] discrepancyMatr;


    GaussMethod(double[][] matrA, double[] vectorB){
        n = matrA.length;
        mainElements = new ArrayList<>();
        this.matrA = new double[n][n];
        this.resultMatrA = new double[n][n];
        this.vectorB = new double[n];
        this.copyVectorB = new double[n];
        this.vectorX = new double[n];
        this.discrepancy = new double[n];
        this.queue = new int[2][n];
        for(int i = 0; i < n; i++){
            this.queue[0][i] = i;
            this.queue[1][i] = i;
            this.vectorB[i] = vectorB[i];
            this.copyVectorB[i] = vectorB[i];
            for(int j = 0; j < n; j++){
                this.resultMatrA[i][j] = matrA[i][j];
                this.matrA[i][j] = matrA[i][j];
            }
        }
    }
    private Pair maxElement(int count){
        double max = Math.abs(resultMatrA[count][count]);
        int firstIndex = count;
        int secondIndex = count;
        for(int i = count; i < n; i++){
            for(int j = count; j < n; j++){
                if(abs(resultMatrA[i][j]) > max){
                    max = resultMatrA[i][j];
                    firstIndex = i;
                    secondIndex = j;
                }
            }
        }
        mainElements.add(resultMatrA[firstIndex][secondIndex]);
        return new Pair(firstIndex, secondIndex);
    }
    private void divizionOnElem(int count){
        double elem = resultMatrA[count][count];
        double elem1;
        vectorB[count] /= elem;
        for(int i = count; i < n; i++) {
            resultMatrA[count][i] /= elem;
        }
        for(int i = count + 1; i< n; i++){
            elem1 = resultMatrA[i][count];
            vectorB[i] -= vectorB[count] * elem1;
            for(int j = count; j < n; j++){
                resultMatrA[i][j] -= resultMatrA[count][j]*elem1;
            }
        }
    }


    public void strightStep(){
        for(int count = 0; count < n; count++){
            swapColumnsAndRows(count, maxElement(count));
            divizionOnElem(count);
        }
    }

    private void swapColumnsAndRows(int count, Pair index) {
        double tmp;
        if (count != index.getFirst()) {
            numberOfPermutations++;
            for (int i = count; i < n; i++) {
                tmp = resultMatrA[index.getFirst()][i];
                resultMatrA[index.getFirst()][i] = resultMatrA[count][i];
                resultMatrA[count][i] = tmp;
            }
            double tmp2 = vectorB[index.getFirst()];
            vectorB[index.getFirst()] = vectorB[count];
            vectorB[count] = tmp2;

        }

        if (count != index.getSecond()) {
            numberOfPermutations++;
            for (int j = 0; j < n; j++) {
                tmp = resultMatrA[j][count];
                resultMatrA[j][count] = resultMatrA[j][index.getSecond()];
                resultMatrA[j][index.getSecond()] = tmp;
            }
            int tmp3 = queue[1][index.getSecond()];
            queue[1][index.getSecond()] = queue[1][count];
            queue[1][count] = tmp3;
        }
    }
    public void print(){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                System.out.format("%25s", resultMatrA[i][j] + "    ");
            }
            System.out.print("    " + vectorB[i] + "\n");
        }
    }
    public void createVectorX(){
        double sum;
        vectorX[n-1] = vectorB[n-1];
        for(int i = n - 2; i >= 0; i--){
            sum = 0;
            for(int j = n-1; j > i; j--){
                sum += resultMatrA[i][j] * vectorX[j];
            }
            vectorX[i] = vectorB[i] - sum;
        }
        double tmp[] = new double[n];
        for(int i = 0; i < n; i++) {
            tmp[i] = vectorX[i];
        }
        for(int i = 0; i < n; i++){
            vectorX[queue[1][i]] = tmp[i];
        }
    }

    public void printX(){
        for(int i = 0; i < n; i++){
            System.out.println(vectorX[i]);
        }
    }

    public void createDiscrepancy(){
        double[] res = new double[n];
        for(int i = 0; i < n; i++){
            res[i] = 0;
            for(int j = 0; j < n; j++){
                res[i] += matrA[i][j]*vectorX[j];
            }
            discrepancy[i] = copyVectorB[i] - res[i];
            System.out.println(discrepancy[i]);
        }
    }

    public void createDetA(){
        detA = Math.pow(-1, numberOfPermutations);
        for(Double item: mainElements){
            detA *= item;
        }
        System.out.println("\ndet A = " + detA);
    }

    public double[] getVectorX() {
        return vectorX;
    }

    public void createDiscrepancyMatr(double inverseMatr[][]){
        discrepancyMatr = new double[n][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                discrepancyMatr[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    discrepancyMatr[i][j] += matrA[i][k] * inverseMatr[k][j];
                }
                if(i == j){
                    discrepancyMatr[i][j] -= 1;
                }
            }
        }
    }

    public void printDiscrepancyMatr(){
        System.out.println("\n Discrepancy matr:");
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                System.out.format("%25s", discrepancyMatr[i][j] + "    ");
            }
            System.out.println();
        }
    }

}

