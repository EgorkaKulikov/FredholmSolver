using System;
using System.Collections.Generic;

namespace FredholmSolver
{
  public class FunctionalConstructor
  {
    /// <summary>
    /// Вычисление значения аппроксимационного функционала.
    /// 
    /// Значение определяется типом и номером выбранного функционала, 
    /// функцией, к которой он примененён, а также основной и дополнительной сеткой.
    /// </summary>
    public double FunctionalValue(
      int funcNumber,
      Func<double, double> f,
      Dictionary<int, double> grid,
      Dictionary<int, double> additionalGrid
      )
    {
      int j = funcNumber;
      
      Func<double, double> Derivative(Func<double, double> func, double h = 1e-5) => 
        x => (func(x + h) - func(x - h)) / (2 * h);
      
      Func<double, double> fDer = Derivative(f);
      Func<double, double> fDer2 = Derivative(fDer);

      switch (Configuration.ApproxType)
      {
        case ApproximationType.ShoenbergMarsden:
          return f(additionalGrid[j]);

        case ApproximationType.Averaging:
          if (j == -2 || j == Configuration.GridPoints - 1)
            return f(additionalGrid[j]);

          double a = -1.0 / 8;
          double b = 5.0 / 4;
          double c = -1.0 / 8;
          return a * f(additionalGrid[j - 1]) + b * f(additionalGrid[j]) + c * f(additionalGrid[j + 1]);
        
        case ApproximationType.DeBoorFix0:
          return f(grid[j]) + (1.0 / 2 * (grid[j + 1] + grid[j + 2]) - grid[j]) * fDer(grid[j])
                            + 1.0 / 2 * (grid[j + 1] + grid[j]) * (grid[j + 2] - grid[j]) * fDer2(grid[j]);
        
        case ApproximationType.DeBoorFix1:
          return f(grid[j + 1]) + 1.0 / 2 * (grid[j + 2] - grid[j + 1]) * fDer(grid[j + 1]);
        
        case ApproximationType.DeBoorFix2:
          return f(grid[j + 2]) - 1.0 / 2 * (grid[j + 2] - grid[j + 1]) * fDer(grid[j + 2]);
        
        case ApproximationType.Proectional:
          if (j == -2)
            return f(grid[0]);
          
          if (j == Configuration.GridPoints - 1)
            return f(grid[Configuration.GridPoints]);
          
          if (j == -1)
            return -1.0 / 2 * f(grid[0]) + 2 * f(0.5 * (grid[0] + grid[1])) - 1.0 / 2 * f(grid[1]);
          
          if (j == Configuration.GridPoints - 2)
            return -1.0 / 2 * f(grid[Configuration.GridPoints - 1])
              + 2 * f(0.5 * (grid[Configuration.GridPoints - 1] + grid[Configuration.GridPoints])) 
              -1.0 / 2 * f(grid[Configuration.GridPoints]);


          return 1.0 / 14 * f(grid[j]) 
                 - 2.0 / 7 * f(0.5 * (grid[j] + grid[j + 1]))
                 + 10.0 / 7 * f(0.5 * (grid[j + 1] + grid[j + 2]))
                 - 2.0 / 7 * f(0.5 * (grid[j + 2] + grid[j + 3]))
                 + 1.0 / 14 * f(grid[j + 3]);

        case ApproximationType.AveragingTrigonometric:
          if (j == -2 || j == Configuration.GridPoints - 1)
          {
            return f(additionalGrid[j]);
          }

          double S1 = 1.0 * (Math.Cos(grid[j + 1]) - Math.Cos(grid[j + 2])) / Math.Sin(grid[j + 2] - grid[j + 1]);
          double S2 = 1.0 * (Math.Sin(grid[j + 2]) - Math.Sin(grid[j + 1])) / Math.Sin(grid[j + 2] - grid[j + 1]);

          double zn = (Math.Sin(additionalGrid[j + 1]) - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j]) - Math.Cos(additionalGrid[j - 1]))
            - (Math.Sin(additionalGrid[j]) - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j + 1]) - Math.Cos(additionalGrid[j - 1]));

          double ch_b = (Math.Sin(additionalGrid[j + 1]) - Math.Sin(additionalGrid[j - 1])) * (S2 - Math.Cos(additionalGrid[j - 1]))
            - (S1 - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j + 1]) - Math.Cos(additionalGrid[j - 1]));

          double ch_c = (S1 - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j]) - Math.Cos(additionalGrid[j - 1]))
            - (Math.Sin(additionalGrid[j]) - Math.Sin(additionalGrid[j - 1])) * (S2 - Math.Cos(additionalGrid[j - 1]));

          b = 1.0 * ch_b / zn;
          c = 1.0 * ch_c / zn;
          a = 1 - b - c;
          return a * f(additionalGrid[j - 1]) + b * f(additionalGrid[j]) + c * f(additionalGrid[j + 1]);

        case ApproximationType.AveragingHyperbolic:
          if (j == -2 || j == Configuration.GridPoints - 1)
          {
            return f(additionalGrid[j]);
          }

          S1 = 1.0 * (Math.Cosh(grid[j + 2]) - Math.Cosh(grid[j + 1])) / Math.Sinh(grid[j + 2] - grid[j + 1]);
          S2 = 1.0 * (Math.Sinh(grid[j + 2]) - Math.Sinh(grid[j + 1])) / Math.Sinh(grid[j + 2] - grid[j + 1]);

          zn = (Math.Sinh(additionalGrid[j + 1]) - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j]) - Math.Cosh(additionalGrid[j - 1]))
            - (Math.Sinh(additionalGrid[j]) - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j + 1]) - Math.Cosh(additionalGrid[j - 1]));

          ch_b = (Math.Sinh(additionalGrid[j + 1]) - Math.Sinh(additionalGrid[j - 1])) * (S2 - Math.Cosh(additionalGrid[j - 1]))
            - (S1 - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j + 1]) - Math.Cosh(additionalGrid[j - 1]));

          ch_c = (S1 - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j]) - Math.Cosh(additionalGrid[j - 1]))
            - (Math.Sinh(additionalGrid[j]) - Math.Sinh(additionalGrid[j - 1])) * (S2 - Math.Cosh(additionalGrid[j - 1]));

          b = 1.0 * ch_b / zn;
          c = 1.0 * ch_c / zn;
          a = 1 - b - c;
          return a * f(additionalGrid[j - 1]) + b * f(additionalGrid[j]) + c * f(additionalGrid[j + 1]);
        
        case ApproximationType.ProectionalHyperbolic:
        case ApproximationType.ProectionalTrigonometric:
          if (j == -2)
            return f(grid[0]);
          
          if (j == Configuration.GridPoints - 1)
            return f(grid[Configuration.GridPoints]);
          
          SplineConstructor sc = new SplineConstructor();
          var A = sc.SplineValueAtPoint(0, grid, grid[1]);
          var B = sc.SplineValueAtPoint(0, grid, 0.5 * (grid[0] + grid[1]));
          var C = sc.SplineValueAtPoint(0, grid, 0.5 * (grid[1] + grid[2]));
          var D = sc.SplineValueAtPoint(0, grid, grid[2]);
          var E = sc.SplineValueAtPoint(0, grid, 0.5 * (grid[2] + grid[3]));
          
          double K = -1.0 / ( A * C * D - A * A * E - B * D * D);
          double K1 = A * E * K;
          double K2 = -A * D * K;
          double K3 = B * D * K;
          
          if (j == -1)
            return K1 * f(grid[0]) + K2 * f(0.5 * (grid[0] + grid[1])) + K3 * f(grid[1]);
          
          if (j == Configuration.GridPoints - 2)
            return K1 * f(grid[Configuration.GridPoints - 1])
                   + K2 * f(0.5 * (grid[Configuration.GridPoints - 1] + grid[Configuration.GridPoints])) 
                   + K3 * f(grid[Configuration.GridPoints]);
          
          K = 1.0 / (C * C * D - A * C * E - B * D * E - D * E * E);
          K1 = E * E * K;
          K2 = - D * E * K;
          K3 = (C * D - A * E) * K;


          return K1 * f(grid[j]) 
                 + K2 * f(0.5 * (grid[j] + grid[j + 1]))
                 + K3 * f(0.5 * (grid[j + 1] + grid[j + 2]))
                 + K2 * f(0.5 * (grid[j + 2] + grid[j + 3]))
                 + K1 * f(grid[j + 3]);

        default:
          throw new NotSupportedException($"Functional type {Configuration.ApproxType.ToString()} is not supported");
      }
    }
  }
}
