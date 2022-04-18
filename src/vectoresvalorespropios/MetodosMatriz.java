package vectoresvalorespropios;

import java.text.DecimalFormat;

public class MetodosMatriz implements Cloneable {

    DecimalFormat df = new DecimalFormat("0.000");
    public int n;      //dimensión
    public double[][] x;
    double[][] matrizIm;
    double[][] matrizVec;

    public MetodosMatriz(int n) {
        this.n = n;
        x = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i][j] = 0.0;
            }
        }
    }

    public MetodosMatriz(double[][] x) {
        this.x = x;
        n = x.length;
    }

    public Object clone() {
        MetodosMatriz obj = null;
        try {
            obj = (MetodosMatriz) super.clone();
        } catch (CloneNotSupportedException ex) {
            System.out.println(" no se puede duplicar");
        }
//aquí está la clave  para clonar la matriz bidimensional
        obj.x = (double[][]) obj.x.clone();
        for (int i = 0; i < obj.x.length; i++) {
            obj.x[i] = (double[]) obj.x[i].clone();
        }
        return obj;
    }

    double traza() {
        double tr = 0.0;
        for (int i = 0; i < n; i++) {
            tr += x[i][i];
        }
        return tr;
    }
//producto de dos matrices

    static MetodosMatriz producto(MetodosMatriz a, MetodosMatriz b) {
        MetodosMatriz resultado = new MetodosMatriz(a.n);
        for (int i = 0; i < a.n; i++) {
            for (int j = 0; j < a.n; j++) {
                for (int k = 0; k < a.n; k++) {
                    resultado.x[i][j] += a.x[i][k] * b.x[k][j];
                }
            }
        }
        return resultado;
    }
//matriz inversa

    static MetodosMatriz inversa(MetodosMatriz d) {
        int n = d.n;  //dimensión de la matriz
        MetodosMatriz a = (MetodosMatriz) d.clone();
        MetodosMatriz b = new MetodosMatriz(n);   //matriz de los términos independientes
        MetodosMatriz c = new MetodosMatriz(n);   //matriz de las incógnitas
//matriz unidad
        for (int i = 0; i < n; i++) {
            b.x[i][i] = 1.0;
        }
//transformación de la matriz y de los términos independientes
        for (int k = 0; k < n - 1; k++) {
            for (int i = k + 1; i < n; i++) {
//términos independientes
                for (int s = 0; s < n; s++) {
                    b.x[i][s] -= a.x[i][k] * b.x[k][s] / a.x[k][k];
                }
//elementos de la matriz
                for (int j = k + 1; j < n; j++) {
                    a.x[i][j] -= a.x[i][k] * a.x[k][j] / a.x[k][k];
                }
            }
        }
//cálculo de las incógnitas, elementos de la matriz inversa
        for (int s = 0; s < n; s++) {
            c.x[n - 1][s] = b.x[n - 1][s] / a.x[n - 1][n - 1];
            for (int i = n - 2; i >= 0; i--) {
                c.x[i][s] = b.x[i][s] / a.x[i][i];
                for (int k = n - 1; k > i; k--) {
                    c.x[i][s] -= a.x[i][k] * c.x[k][s] / a.x[i][i];
                }
            }
        }
        return c;
    }
//matriz traspuesta

    static MetodosMatriz traspuesta(MetodosMatriz a) {
        int n = a.n;    //dimensión
        MetodosMatriz d = new MetodosMatriz(a.n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                d.x[i][j] = a.x[j][i];
            }
        }
        return d;
    }

//Valores propios 
   public MetodosMatriz valoresPropios(double[] valores, int maxIter) throws ValoresExcepcion {
        final double CERO = 1e-8;
        double maximo, tolerancia, sumsq;
        double x, y, z, c, s;
        int contador = 0;
        int i, j, k, l;
        MetodosMatriz a = (MetodosMatriz) clone();      //matriz copia
        MetodosMatriz p = new MetodosMatriz(n);
        do {
            k = 0;
            l = 1;
            maximo = Math.abs(a.x[k][1]);
            for (i = 0; i < n - 1; i++) {
                for (j = i + 1; j < n; j++) {
                    if (Math.abs(a.x[i][j]) > maximo) {
                        k = i;
                        l = j;
                        maximo = Math.abs(a.x[i][j]);
                    }
                }
            }
            sumsq = 0.0;
            for (i = 0; i < n; i++) {
                sumsq += a.x[i][i] * a.x[i][i];
            }
            tolerancia = 0.0001 * Math.sqrt(sumsq) / n;
            if (maximo < tolerancia) {
                break;
            }
//calcula la matriz ortogonal de p
//inicialmente es la matriz unidad
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    p.x[i][j] = 0.0;
                }
            }
            for (i = 0; i < n; i++) {
                p.x[i][i] = 1.0;
            }
            y = a.x[k][k] - a.x[l][l];
            if (Math.abs(y) < CERO) {
                c = s = Math.sin(Math.PI / 4);
            } else {
                x = 2 * a.x[k][l];
                z = Math.sqrt(x * x + y * y);
                c = Math.sqrt((z + y) / (2 * z));
                s = signo(x / y) * Math.sqrt((z - y) / (2 * z));
            }
            p.x[k][k] = c;
            p.x[l][l] = c;
            p.x[k][l] = s;
            p.x[l][k] = -s;
            a = MetodosMatriz.producto(p, MetodosMatriz.producto(a, MetodosMatriz.traspuesta(p)));
            MatrizP.Rvalores.append("--------------------------------------------------------------------------");
            MatrizP.Rvalores.append(a.mostrarDatos());
            contador++;
        } while (contador < maxIter);
        if (contador == maxIter) {
            throw new ValoresExcepcion("No se han podido calcular los valores propios");
        }
//valores propios
        for (i = 0; i < n; i++) {
            valores[i] = (double) Math.round(a.x[i][i] * 1000) / 1000;
        }

        return a;
    }

    public MetodosMatriz vectoresPropios(double[] valores, int maxIter) throws ValoresExcepcion {
        final double CERO = 1e-8;
        double maximo, tolerancia, sumsq;
        double x, y, z, c, s;
        int contador = 0;
        int i, j, k, l;
        MetodosMatriz a = (MetodosMatriz) clone();      //matriz copia
        MetodosMatriz p = new MetodosMatriz(n);
        MetodosMatriz q = new MetodosMatriz(n);
//matriz unidad
        for (i = 0; i < n; i++) {
            q.x[i][i] = 1.0;
        }
        do {
            k = 0;
            l = 1;
            maximo = Math.abs(a.x[k][1]);
            for (i = 0; i < n - 1; i++) {
                for (j = i + 1; j < n; j++) {
                    if (Math.abs(a.x[i][j]) > maximo) {
                        k = i;
                        l = j;
                        maximo = Math.abs(a.x[i][j]);
                    }
                }
            }
            sumsq = 0.0;
            for (i = 0; i < n; i++) {
                sumsq += a.x[i][i] * a.x[i][i];
            }
            tolerancia = 0.0001 * Math.sqrt(sumsq) / n;
            if (maximo < tolerancia) {
                break;
            }
//calcula la matriz ortogonal de p
//inicialmente es la matriz unidad
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    p.x[i][j] = 0.0;
                }
            }
            for (i = 0; i < n; i++) {
                p.x[i][i] = 1.0;
            }
            
            y = a.x[k][k] - a.x[l][l];
            if (Math.abs(y) < CERO) {
                c = s = Math.sin(Math.PI / 4);
            } else {
                x = 2 * a.x[k][l];
                z = Math.sqrt(x * x + y * y);
                c = Math.sqrt((z + y) / (2 * z));
                s = signo(x / y) * Math.sqrt((z - y) / (2 * z));
            }

            p.x[k][k] = c;
            p.x[l][l] = c;
            p.x[k][l] = s;
            p.x[l][k] = -s;
            a = MetodosMatriz.producto(p, MetodosMatriz.producto(a, MetodosMatriz.traspuesta(p)));
            q = MetodosMatriz.producto(q, MetodosMatriz.traspuesta(p));
            /////////////////////////////////
            matrizIm = new double[n][n];

            for (int b = 0; b < matrizIm.length; b++) {
                for (int m = 0; m < matrizIm.length; m++) {
                    matrizIm[b][m] = q.x[b][m];
                }
            }

            MatrizP.Rvectores.append("--------------------------------------------------------------------------\n");
            for (int b = 0; b < matrizIm.length; b++) {
                for (int m = 0; m < matrizIm.length; m++) {
                    MatrizP.Rvectores.append(df.format(matrizIm[b][m]) + "\t");
                }
                MatrizP.Rvectores.append("\n");
            }
            ///////////////////////
            contador++;
        } while (contador < maxIter);

        if (contador == maxIter) {
            throw new ValoresExcepcion("No se han podido calcular los valores propios");
        }
//valores propios
        for (i = 0; i < n; i++) {
            valores[i] = (double) Math.round(a.x[i][i] * 1000) / 1000;
        }
//vectores propios
        return q;
    }

    public double rd(double valorInicial, int numeroDecimales) {
        double parteEntera, resultado;
        resultado = valorInicial;
        parteEntera = Math.floor(resultado);
        resultado = (resultado - parteEntera) * Math.pow(10, numeroDecimales);
        resultado = Math.round(resultado);
        resultado = (resultado / Math.pow(10, numeroDecimales)) + parteEntera;
        return resultado;
    }

    private int signo(double x) {
        return (x > 0 ? 1 : -1);
    }

    public String mostrarDatos() {
        String texto = "\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                texto += ((double) Math.round(1000 * x[i][j]) / 1000) + "\t";
            }
            texto += "\n";
        }
        texto += "\n";
        return texto;
    }
}

class ValoresExcepcion extends Exception {

    public ValoresExcepcion() {
        super();
    }

    public ValoresExcepcion(String s) {
        super(s);
    }
}
