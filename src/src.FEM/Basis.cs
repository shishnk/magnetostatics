namespace Magnetostatics.src.FEM;

public interface IBasis2D
{
    int Size { get; }

    double GetPsi(int number, Point2D point);

    double GetDPsi(int number, int varNumber, Point2D point);
}

public interface IBasis1D
{
    int Size { get; }

    double GetPsi(int number, double point, double h = 1.0);

    double GetDPsi(int number, double point, double h = 1.0);

    double GetDdPsi(int number, double point, double h = 1.0);
}

public class HermiteBasis : IBasis1D
{
    public int Size => 4;

    public double GetPsi(int number, double point, double h = 1.0)
        => number switch
        {
            0 => 1.0 - 3.0 * point * point + 2.0 * point * point * point, // $1 - 3\xi^2 + 2\xi^3$
            1 => h * (point - 2.0 * point * point + point * point * point), // $h_i \cdot (\xi - 2\xi^2 + \xi^3)$
            2 => 3.0 * point * point - 2.0 * point * point * point, // $3\xi^2 - 2\xi^3$
            3 => h * (-point * point + point * point * point), // $h_i \cdot (-\xi^2 + \xi^3)$
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
        };

    public double GetDPsi(int number, double point, double h = 1.0)
        => number switch
        {
            0 => -6.0 * (point - point * point) / h, // $\dfrac{-6 \cdot (\xi - \xi^2)}{h_i}$
            1 => 1.0 - 4 * point + 3 * point * point, // $1 - 4\xi + 3\xi^2$
            2 => 6.0 * (point - point * point) / h, // $\dfrac{6 \cdot (\xi - \xi^2)}{h_i}$
            3 => -2.0 * point + 3 * point * point, // $-2\xi + 3\xi^2$
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
        };

    public double GetDdPsi(int number, double point, double h = 1.0)
        => number switch
        {
            0 => -6.0 * (1.0 - 2.0 * point) / (h * h), // $\dfrac{-6 \cdot (1 - 2\xi)}{h_i^2}$,
            1 => (-4.0 + 6.0 * point) / h, // $\dfrac{-4 + 6\xi}{h_i}$
            2 => 6.0 * (1.0 - 2.0 * point) / (h * h), // $\dfrac{6 \cdot (1 - 2\xi)}{h_i^2}$
            3 => (-2.0 + 6.0 * point) / h, // $\dfrac{-2 + 6\xi}{h_i}$
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
        };
}

public readonly record struct LinearBasis : IBasis2D
{
    public int Size => 4;

    public double GetPsi(int number, Point2D point)
        => number switch
        {
            0 => (1.0 - point.X) * (1.0 - point.Y),
            1 => point.X * (1.0 - point.Y),
            2 => (1.0 - point.X) * point.Y,
            3 => point.X * point.Y,
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
        };

    public double GetDPsi(int number, int varNumber, Point2D point)
        => varNumber switch
        {
            0 => number switch
            {
                0 => point.Y - 1.0,
                1 => 1.0 - point.Y,
                2 => -point.Y,
                3 => point.Y,
                _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
            },
            1 => number switch
            {
                0 => point.X - 1.0,
                1 => -point.X,
                2 => 1.0 - point.X,
                3 => point.X,
                _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
            },
            _ => throw new ArgumentOutOfRangeException(nameof(varNumber), varNumber, "Not expected var number!")
        };
}