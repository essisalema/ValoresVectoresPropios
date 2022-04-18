package vectoresvalorespropios;

import java.text.DecimalFormat;
import javax.swing.JOptionPane;
import javax.swing.table.DefaultTableModel;

public class MatrizP extends javax.swing.JFrame {

    static int f = 0, c = 0;
    DefaultTableModel modelo;
    double[][] matriz;
    double[][] matrizVectores;
    double[][] matrizIm;
    double[][] diagonal;
    DecimalFormat df = new DecimalFormat("0.000");

    /**
     * Creates new form matriz
     */
    public MatrizP() {
        initComponents();
        modelo = new DefaultTableModel(f, c);
        datosMatriz.setModel(modelo);
        this.setLocationRelativeTo(null);
        this.setResizable(false);
        this.setTitle("Valores y Vectores Propios");

    }

    public static void dimension() {
        try {
            f = Integer.parseInt(JOptionPane.showInputDialog(null, "Ingrese el numero de filas de la matriz"));
            c = Integer.parseInt(JOptionPane.showInputDialog(null, "Ingrese el numero de columnas de la matriz"));
        } catch (NumberFormatException e) {
            JOptionPane.showMessageDialog(null, "ERROR: " + e.getMessage());
            System.exit(0);
        }
    }

    public void obtenerDatos() {
        matriz = new double[f][c];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz.length; j++) {
                matriz[i][j] = Double.parseDouble(datosMatriz.getValueAt(i, j).toString());
            }
        }
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane1 = new javax.swing.JScrollPane();
        datosMatriz = new javax.swing.JTable();
        btnValores = new javax.swing.JButton();
        btnVectores = new javax.swing.JButton();
        jScrollPane2 = new javax.swing.JScrollPane();
        Rvectores = new javax.swing.JTextArea();
        btnDimension = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        jScrollPane3 = new javax.swing.JScrollPane();
        Rvalores = new javax.swing.JTextArea();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        datosMatriz.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {},
                {},
                {},
                {},
                {}
            },
            new String [] {

            }
        ));
        jScrollPane1.setViewportView(datosMatriz);

        btnValores.setText("Valores Propios");
        btnValores.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnValoresActionPerformed(evt);
            }
        });

        btnVectores.setText("Vectores Propios");
        btnVectores.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnVectoresActionPerformed(evt);
            }
        });

        Rvectores.setColumns(20);
        Rvectores.setRows(5);
        jScrollPane2.setViewportView(Rvectores);

        btnDimension.setText("Nueva Dimensión");
        btnDimension.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnDimensionActionPerformed(evt);
            }
        });

        jLabel1.setFont(new java.awt.Font("Bookman Old Style", 0, 14)); // NOI18N
        jLabel1.setText("Ingrese los datos");

        Rvalores.setColumns(20);
        Rvalores.setRows(5);
        jScrollPane3.setViewportView(Rvalores);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 375, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 375, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jLabel1, javax.swing.GroupLayout.PREFERRED_SIZE, 144, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(66, 66, 66)
                        .addComponent(btnValores)
                        .addGap(48, 48, 48)
                        .addComponent(btnVectores))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jScrollPane3, javax.swing.GroupLayout.PREFERRED_SIZE, 375, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(137, 137, 137)
                        .addComponent(btnDimension)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 113, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnValores)
                    .addComponent(btnVectores))
                .addGap(18, 18, 18)
                .addComponent(btnDimension)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane3, javax.swing.GroupLayout.PREFERRED_SIZE, 133, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 267, Short.MAX_VALUE)
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void btnValoresActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnValoresActionPerformed
        try {
            obtenerDatos();
            Rvalores.setText("");
            MetodosMatriz mm = new MetodosMatriz(matriz);
            MetodosMatriz valoresP = new MetodosMatriz(mm.n);
            double[] valores = new double[mm.n];
            try {
                valoresP = mm.valoresPropios(valores, 100);

            } catch (ValoresExcepcion ex) {
                JOptionPane.showMessageDialog(null, "Al calcular los valores propios se ha producido una excepción\n "
                        + ex.getClass() + " con el mensaje " + ex.getMessage());
            }
            Rvalores.append("--------------------------------------------------------------------------\n");
            Rvalores.append("Valores Propios\n");
            String lambda = "\u03BB";
            for (int i = 0; i < mm.n; i++) {            //valores propios
                Rvalores.append(lambda + (i + 1) + ": " + valores[i] + "\n");
            }
        } catch (NullPointerException e) {
            JOptionPane.showMessageDialog(null, "Ingrese los datos");
        }
    }//GEN-LAST:event_btnValoresActionPerformed

    private void btnVectoresActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnVectoresActionPerformed
        try {
            Rvectores.setText("");
            obtenerDatos();
            MetodosMatriz mm = new MetodosMatriz(matriz);
            MetodosMatriz vectores = new MetodosMatriz(mm.n);
            double[] valores = new double[mm.n];
            try {
                vectores = mm.vectoresPropios(valores, 100);
            } catch (ValoresExcepcion ex) {
                JOptionPane.showMessageDialog(null, "Al calcular los valores propios se ha producido una excepción\n "
                        + ex.getClass() + " con el mensaje " + ex.getMessage());
            }

            ////////////////////////////////////////
            MetodosMatriz a = (MetodosMatriz) vectores.clone();
            matrizVectores = new double[f][c];
            diagonal = new double[f][c];
            ///Matriz resultado
            for (int i = 0; i < matrizVectores.length; i++) {
                for (int j = 0; j < matrizVectores.length; j++) {
                    matrizVectores[i][j] = a.x[i][j];
                }
            }
            //Matriz diagonal valores propios        
            for (int i = 0; i < diagonal.length; i++) {
                diagonal[i][i] = valores[i];
            }
            MatrizP.Rvectores.append("--------------------------------------------------------------------------\n");
            Rvectores.append("\nVectores Propios: " + vectores.mostrarDatos());
            comprobacionVectores();
        } catch (NullPointerException e) {
            JOptionPane.showMessageDialog(null, "Ingrese los datos");
        }
    }//GEN-LAST:event_btnVectoresActionPerformed

    public void comprobacionVectores() {
        MetodosMatriz matRes = new MetodosMatriz(matrizVectores);
        MetodosMatriz matDiag = new MetodosMatriz(diagonal);
        matRes = MetodosMatriz.producto(matRes, MetodosMatriz.producto(matDiag, MetodosMatriz.inversa(matRes)));
        Rvectores.append("--------------------------------------------------------------------------\n");
        Rvectores.append("Comprobacion Matriz Inicio: " + matRes.mostrarDatos());
    }
    private void btnDimensionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnDimensionActionPerformed
        // TODO add your handling code here:
        dimension();
        if (f >= 2 && c >= 2) {
            if (f == c) {
                new MatrizP().setVisible(true);
                dispose();
            } else {
                JOptionPane.showMessageDialog(null, "No es una matriz simetrica");
                System.exit(0);
            }
        } else {
            JOptionPane.showMessageDialog(null, "Ha ingresado una dimensión invalida");
            System.exit(0);
        }
    }//GEN-LAST:event_btnDimensionActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;

                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(MatrizP.class
                    .getName()).log(java.util.logging.Level.SEVERE, null, ex);

        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(MatrizP.class
                    .getName()).log(java.util.logging.Level.SEVERE, null, ex);

        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(MatrizP.class
                    .getName()).log(java.util.logging.Level.SEVERE, null, ex);

        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(MatrizP.class
                    .getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
        dimension();
        if (f >= 2 && c >= 2) {
            if (f == c) {
                java.awt.EventQueue.invokeLater(new Runnable() {
                    public void run() {
                        new MatrizP().setVisible(true);
                    }
                });
            } else {
                JOptionPane.showMessageDialog(null, "No es una matriz simetrica");
                System.exit(0);
            }
        } else {
            JOptionPane.showMessageDialog(null, "Ha ingresado una dimensión invalida");
            System.exit(0);
        }

    }


    // Variables declaration - do not modify//GEN-BEGIN:variables
    public static javax.swing.JTextArea Rvalores;
    public static javax.swing.JTextArea Rvectores;
    private javax.swing.JButton btnDimension;
    private javax.swing.JButton btnValores;
    private javax.swing.JButton btnVectores;
    private javax.swing.JTable datosMatriz;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JScrollPane jScrollPane3;
    // End of variables declaration//GEN-END:variables
}