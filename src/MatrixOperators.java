import java.util.HashSet;
import java.util.Set;

public class MatrixOperators {
    public MatrixOperators() {
    }

    public double determinant(double[][] matrix) {
        if (matrix.length == 0) {
            return 1;
        } else if (matrix.length != matrix[0].length) {
            System.out.println("Determinants can only be found for square matrices.");
            return 0;
        }

        double determinant = 0;
        for (int j = 0; j < matrix[0].length; j++) {
            determinant += matrix[0][j] * cofactor(matrix, 0, j);
        }
        return determinant;
    }

    public double[][] matrix_of_CoFtr(double[][] matrix) {
        double[][] m_o_cfs = new double[matrix[0].length][matrix.length];

        for (int i = 0; i < matrix[0].length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                m_o_cfs[i][j] = cofactor(matrix, i, j);
            }
        }
        return m_o_cfs;
    }

    public double[][] adjoint(double[][] matrix) {
        return transpose(matrix_of_CoFtr(matrix));
    }

    public double[][] inverse(double[][] matrix) {
        return scmultiply(1/determinant(matrix), adjoint(matrix));
    }

    public double[][] scmultiply(double d, double[][] matrix) {
        double[][] result = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix[0].length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                result[i][j] = matrix[i][j] * d;
            }
        }
        return result;
    }

    private double[][] submatrix(double[][] matrix, Set<Integer> rows_removed, Set<Integer> cols_removed) {
        double[][] m = new double[matrix.length - rows_removed.size()][matrix[0].length - cols_removed.size()];

        int skip_row = 0;
        int skip_col;
        for (int i = 0; i < m.length; i++) {
            if (rows_removed.contains(i + skip_row)) {
                skip_row ++;
            }
            skip_col = 0;
            for (int j = 0; j < m[0].length; j++) {
                if (cols_removed.contains(j + skip_col)) {
                    skip_col ++;
                }
                m[i][j] = matrix[i + skip_row][j + skip_col];
            }
        }
        return m;
    }

    public void print_matrix(double[][] matrix) {
        // finding longest number
        int[] longests = new int[matrix[0].length];
        for (double[] row : matrix) {
            int col = 0;
            for (double n : row) {
                if ((int) n == n) {
                    if (String.valueOf((int) n).length() > longests[col])
                        longests[col] = String.valueOf((int) n).length();
                } else {
                    if (String.valueOf(n).length() > longests[col])
                        longests[col] = String.valueOf(n).length();
                }
                col ++;
            }
        }
        for (double[] row : matrix) {
            System.out.print("[");
            int col = 0;
            for (double n : row) {
                //space
                System.out.print(" ");
                // number
                if ((int) n == n)
                    System.out.print((int) n);
                else
                    System.out.print(n);

                // spaces
                int lofn;
                if ((int) n == n) {
                    lofn = String.valueOf((int) n).length();
                } else {
                    lofn = String.valueOf(n).length();
                }
                // System.out.println("n = " + n + ", lofn = " + lofn);
                // System.out.println("Longest - lofn = " + (longests[col] - lofn));

                for (int i = 0; i < longests[col] - lofn; i ++) {
                    System.out.print(" ");
                }
                System.out.print(" ");
                col ++;
            }
            System.out.println("]");
        }
    }

    public double[][] transpose(double[][] matrix) {
        double[][] transpose = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transpose[j][i] = matrix[i][j];
            }
        }
        return transpose;
    }

    public double trace(double[][] matrix) {
        if (matrix.length != matrix[0].length) {
            System.out.println("Trace can only be found for square matricies");
            return 0;
        }
        double trace = 0;
        for (int i = 0; i < matrix.length; i++) {
            trace += matrix[i][i];
        }
        return trace;
    }

    private double cofactor(double[][] matrix, int i, int j) {
        return minor(matrix, i, j) * Math.pow(-1, (i + j));
    }

    private double minor(double[][] matrix, int i, int j) {
        Set<Integer> rows_removed = new HashSet<>();
        Set<Integer> cols_removed = new HashSet<>();
        rows_removed.add(i);
        cols_removed.add(j);
        return (determinant(submatrix(matrix, rows_removed, cols_removed)));
    }

    public void print_double(double d) {
        if ((int) d == d) {
            System.out.println((int) d);
        } else {
            System.out.println(d);
        }
    }

    public double[][] row_echelon_form(double[][] matrix, boolean printSteps) {
        if (matrix.length == 0 || matrix[0].length == 0) {
            return matrix;
        }

        double[][] m = new double[matrix.length][matrix[0].length];

        for (int r = 0; r < matrix.length; r ++) {
            m[r] = matrix[r].clone();
        }


        // look across the columns for the first one that isn't all 0
        int thisRow = 0;
        int thisCol = 0;

        while (thisRow < m.length && thisCol < m[0].length) {
            int[] nextSpot = findNextNonzeroSpot(m, thisRow, thisCol);

            if (nextSpot[0] == -1) {  // everything else is 0
                return m;
            }

            thisCol = nextSpot[1];

            switch_rows(m, thisRow, nextSpot[0]);

            if (thisRow != nextSpot[0] && printSteps) {
                System.out.println("Switch rows " + thisRow + " and " + nextSpot[0]);
                print_matrix(m);
            }

            double multiplier = 1d / m[thisRow][nextSpot[1]];
            multiply_a_row(m, thisRow, multiplier);  // create a leading 1

            if (multiplier != 1 && printSteps) {
                System.out.println("Multiply row " + thisRow + " by " + multiplier);
                print_matrix(m);
            }

            for (int r = thisRow + 1; r < m.length; r++) {
                if (m[r][thisCol] != 0) {
                    if (printSteps) {
                        System.out.println("Add " + (-1 * m[r][thisCol]) + " times row " + thisRow + " to row " + r);
                    }

                    add_rows(m, r, thisRow, -1 * m[r][thisCol]); // make thisCol of all the other rows 0

                    if (printSteps) {
                        print_matrix(m);
                    }
                }
            }

            thisRow++;
            thisCol++;
        }

        if (printSteps) {
            System.out.println("One row echelon form matrix coming right up: ");
        }
        return m;
    }

    public double[][] reduced_row_echelon_form(double[][] matrix, boolean printSteps) {
        double[][] m = row_echelon_form(matrix, printSteps);

        if (m.length < 2 || m[0].length < 2) {
            if (printSteps) {
                System.out.println("Matrix was too small: returning REF as RREF.");
            }
            return m;
        }

        int col;
        for (int row = 1; row < m.length; row ++) {  // move down the rows
            col = num_leading_zeros(m[row]);
            if (col == m[0].length) {
                if (printSteps) {
                    System.out.println("Row of 0s found. RREF has been achieved.");
                }
                break;
            }
            for (int r = 0; r < row; r ++) {
                if (printSteps) {
                    System.out.println("Add " + (-1 * m[r][col]) + " times row " + row + " to row " + r);
                }

                add_rows(m, r, row, -1 * m[r][col]);

                if (printSteps) {
                    print_matrix(m);
                }
            }
        }

        if (printSteps) {
            System.out.println("One RREF matrix coming right up ...");
        }
        return m;
    }

    private int[] findNextNonzeroSpot(double[][] matrix, int startingRow, int startingCol) {  // returns [row, col] of the first row that doesn't have a 0 in the first nonzero column of matrix

        for (int c = startingCol; c < matrix[0].length; c ++) {
            for (int r = startingRow; r < matrix.length; r ++) {
                if (matrix[r][c] != 0) {
                    return new int[] {r, c};
                }
            }
        }
        return new int[] {-1, -1};
    }

    private double[][] add_rows(double[][] matrix, int targetRow, int rowToAdd, double rowToAddMultiplier) {
        for (int i = 0; i < matrix[targetRow].length; i ++) {
            matrix[targetRow][i] += rowToAddMultiplier * matrix[rowToAdd][i];
        }
        return matrix;
    }

    private double[][] multiply_a_row(double[][] matrix, int targetRow, double multiplier) {
        for (int i = 0; i < matrix[targetRow].length; i ++) {
            matrix[targetRow][i] *= multiplier;
        }
        return matrix;
    }

    private double[][] switch_rows(double[][] matrix, int row1, int row2) {
        double[] tempRow = matrix[row2].clone();
        matrix[row1] = matrix[row2];
        matrix[row2] = tempRow;
        return matrix;
    }

    private int num_leading_zeros(double[] row) {
        for (int i = 0; i < row.length; i ++) {
            if (row[i] != 0) {
                return i;
            }
        }
        return row.length;
    }

}
