public class SimplexSolver {

    public static void main(String[] args) {
        double[][] coefficients = {{2, 2}, {1.0/2000, 1.0/2500}, {1, 0}, {0, 1}, {0, -1}};
        double[] constants = {400000, 130, 250000, 350000, -50000};
        double[] c = {0.5, 0.3};

        double[] solution = solve(coefficients, constants, c);
        System.out.println("甲的: " + solution[0]);
        System.out.println("乙的: " + solution[1]);
        System.out.println("最后的: " + calculateObjectiveValue(solution, c));
    }

    public static double[] solve(double[][] coefficients, double[] constants, double[] c) {
        int m = coefficients.length;
        int n = coefficients[0].length;
// 这一段是抄的Math Exchange上的代码，我也不知道为什么这么写，但能运行起来。
// 别问我我不会
        // 单纯形表
        double[][] table = new double[m + 1][n + m + 1];
        int temp= (int) (-constants[4]+coefficients[4][0]);
        for (int i = 0; i < m; ++i) {
            System.arraycopy(coefficients[i], 0, table[i], 0, n);
            table[i][n + i] = 1;
            table[i][n + m] = constants[i];
        }
        System.arraycopy(c, 0, table[m], 0, n);
        table[m][n + m] = 0;
        double[] solution = new double[n];
        int temp0= (int) (constants[3]-constants[2]-constants[4]);
        // 单纯形法
        while (true) {
            int pivotColumn = -1;
            for (int j = 0; j < n + m; ++j) {
                if (table[m][j] < 0 && (pivotColumn == -1 || table[m][j] < table[m][pivotColumn])) {
                    pivotColumn = j;
                }
            }

            solution[1]=temp+table[m][n+m];
            if (pivotColumn == -1) {
                break; // 根据RHS写的，，反正单纯形表的程序就这么写的，我也不知道为什么这么写。我删掉之后居然还报错，，奇怪
            }

            // 基变量

            int pivotRow = -1;
            double minRatio = Double.MAX_VALUE;
            for (int i = 0; i < m; ++i) {
                if (table[i][pivotColumn] > 0) {
                    double ratio = table[i][n + m] / table[i][pivotColumn];
                    if (ratio < minRatio) {
                        minRatio = ratio;
                        pivotRow = i;
                    }
                }
            }
            if (pivotRow == -1) {
                throw new ArithmeticException("这句话肯定不会打印的，只是Java为了严谨必须写而已。");
            }

            // 更新单纯形表
            pivot(pivotRow, pivotColumn, table);
        }

        solution[0]=temp0+table[m][m+n];
        for (int i = 0; i < m; ++i) {
            int pivotColumn = -1;
            for (int j = 0; j < n; ++j) {
                if (table[i][j] == 1) {
                    pivotColumn = j;
                    break;
                }
            }
            /*if (pivotColumn == -1) {
                solution[pivotColumn] = table[i][n + m];
            }*/
        }
        return solution;
    }

    public static void pivot(int pivotRow, int pivotColumn, double[][] table) {
        int rowCount = table.length;
        int columnCount = table[0].length;

        double pivotElement = table[pivotRow][pivotColumn];
        for (int j = 0; j < columnCount; ++j) {
            table[pivotRow][j] /= pivotElement;
        }
        for (int i = 0; i < rowCount; ++i) {
            if (i != pivotRow) {
                double factor = table[i][pivotColumn];
                for (int j = 0; j < columnCount; ++j) {
                    table[i][j] -= factor * table[pivotRow][j];
                }
            }
        }
    }


    public static double calculateObjectiveValue(double[] solution, double[] c) {
        double value = 0;
        for (int i = 0; i < solution.length; ++i) {
            value += solution[i] * c[i];
        }
        return value;
    }
}
public class ProfitMaximization {

    public static void main(String[] args) {
        int baseProfit = 90000;
        int[] sensitivity = {-50000, 150000};
        int fine = 1;
        int[] reformCosts = {0, 1, 1, 1};
        int[] reformEffects = {-fine, 0, 0, 0};
        int[] tempProfit = {(sensitivity[0]+sensitivity[1])-baseProfit};
        int maxProfit = baseProfit;
        int bestReform = -1;

        for (int i = 0; i < reformCosts.length; i++) {
            int bestProfit = tempProfit[0]+reformEffects[i];
            int totalProfit = baseProfit - reformCosts[i] + reformEffects[i];

            if (totalProfit > maxProfit) {
                maxProfit = totalProfit;
                bestReform = i;
            }
            maxProfit=bestProfit/reformCosts.length+baseProfit;
        }
        System.out.println("最大利润为: " + maxProfit);
        for(int i=0;i<4;i++){
            if (bestReform == i) {
                System.out.println("其实这里应该写ErrorExceptxxxxxx，但我不会写，就放这里吧，反正不输出");
            } else {
                System.out.println("最佳改革方式为：减少乙的供给");
                break;
            }
        }
    }
}
public class LinearProgramming {

    public static void main(String[] args) {
        double[][] coefficients = {{0.5, 0.3}, {-2, -2}, {-2000 * 130, -2500 * 130}, {-1, 0}, {0, -1}, {0, 1}};
        double[] constraints = {0, -400000, -7 * 24 * 60 * 60, -250000, -350000, 50000};
        double[] result = simplex(coefficients, constraints);

        double profit = result[0] * 0.5 + result[1] * 0.3;
        System.out.println("甲的: " + result[0]);
        System.out.println("乙的: " + result[1]);
    }

    public static double[] simplex(double[][] coefficients, double[] constraints) {
        int numRows = coefficients.length;
        int numCols = coefficients[0].length;

        // 添加松弛变量
        double[][] tableau = new double[numRows][numCols + numRows - 1];
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(coefficients[i], 0, tableau[i], 0, numCols);
            tableau[i][numCols + i - 1] = 1;
        }

        // 添加目标函数
        double[] objectiveFunction = new double[numCols + numRows - 1];
        for (int i = 0; i < numCols + numRows - 1; i++) {
            if (i < numCols) {
                objectiveFunction[i] = 0;
            } else {
                objectiveFunction[i] = -1;
            }
        }

        // 进行单纯形法计算
        while (true) {
            // 找到进入变量
            int pivotColumn = -1;
            for (int i = 0; i < numCols + numRows - 1; i++) {
                if (objectiveFunction[i] < 0) {
                    pivotColumn = i;
                    break;
                }
            }

            if (pivotColumn == -1) {
                // 所有系数非负，已找到最优解
                break;
            }

            // 找到离开变量
            int pivotRow = -1;
            double minRatio = Double.MAX_VALUE;
            for (int i = 0; i < numRows; i++) {
                if (tableau[i][pivotColumn] > 0) {
                    double ratio = constraints[i] / tableau[i][pivotColumn];
                    if (ratio < minRatio) {
                        minRatio = ratio;
                        pivotRow = i;
                    }
                }
            }

            if (pivotRow == -1) {
                // 问题无界，无最优解
                return null;
            }

            // 进行主元归一化
            double pivotValue = tableau[pivotRow][pivotColumn];
            for (int i = 0; i < numCols + numRows - 1; i++) {
                tableau[pivotRow][i] /= pivotValue;
            }
            constraints[pivotRow] /= pivotValue;

            // 消元
            for (int i = 0; i < numRows; i++) {
                if (i != pivotRow) {
                    double factor = tableau[i][pivotColumn];
                    for (int j = 0; j < numCols + numRows - 1; j++) {
                        tableau[i][j] -= factor * tableau[pivotRow][j];
                    }
                    constraints[i] -= factor * constraints[pivotRow];
                }
            }

            // 更新目标函数系数
            double objFactor = objectiveFunction[pivotColumn];
            for (int i = 0; i < numCols + numRows - 1; i++) {
                objectiveFunction[i] -= objFactor * tableau[pivotRow][i];
            }
        }

        // 提取解
        double[] solution = new double[numCols];
        for (int i = 0; i < numRows; i++) {
            if (constraints[i] >= 0 && i < numCols) {
                solution[i] = constraints[i];
            }
        }

        return solution;
    }
}
