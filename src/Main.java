import java.util.*;

public class Main {
    public static void main(String[] args) {

        MatrixOperators opr = new MatrixOperators();
        System.out.println("--- Matrix Calculator ---");
        // cols
        Scanner s = new Scanner(System.in);
        int rownum = 0;
        boolean valid_in = false;
        while (!valid_in) {
            System.out.print("Enter the number of rows in your matrix: ");
            try {
                rownum = Integer.parseInt(s.nextLine().strip());
                valid_in = true;
            } catch (NumberFormatException e) {
                System.out.println("The input was invalid. Please try again.");
            }
        }
        double[][] matrix = new double[0][0];
        int colnum;
        // rows
        System.out.println("Enter each entry in your matrix separated by spaces, with 'enter'/'return' when a new row starts:");
        for (int i = 0; i < rownum; i++) {
            valid_in = false;
            while (!valid_in) {
                try {
                    String raw_row_in = s.nextLine();
                    String[] row_in = raw_row_in.split(" ");
                    double[] row = new double[row_in.length];
                    for (int j = 0; j < row_in.length; j++) {
                        row[j] = Double.parseDouble(row_in[j]);
                    }
                    valid_in = true;
                    if (i == 0) {
                        colnum = row_in.length;
                        matrix = new double[rownum][colnum];
                    }
                    matrix[i] = row;
                } catch (NumberFormatException e) {
                    System.out.println("One or more entries were invalid. Please try again:");
                }
            }
        }


        opr.print_matrix(matrix);
        System.out.println("RREF: ");
        opr.print_matrix(opr.reduced_row_echelon_form(matrix, true));

        System.out.print("Determinant: ");
        opr.print_double(opr.determinant(matrix));
        System.out.println("Adjoint:");
        opr.print_matrix(opr.adjoint(matrix));
        System.out.println("Inverse:");
        opr.print_matrix(opr.inverse(matrix));
        System.out.println("Transpose: ");
        opr.print_matrix(opr.transpose(matrix));
        System.out.print("Trace: ");
        opr.print_double(opr.trace(matrix));

    }
}