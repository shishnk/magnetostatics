namespace Magnetostatics.Mesh;

public interface IBaseMesh
{
    IReadOnlyList<Point2D> Points { get; }
    IReadOnlyList<FiniteElement> Elements { get; }
    IReadOnlyList<Area> Areas { get; }
}

public abstract class MeshBuilder
{
    protected abstract int ElementSize { get; }

    public abstract (IReadOnlyList<Point2D>, FiniteElement[]) Build(MeshParameters meshParameters);

    protected (IReadOnlyList<Point2D>, FiniteElement[]) BaseBuild(MeshParameters parameters)
    {
        var result = new
        {
            Points = new Point2D[(parameters.SplitsX.Sum() + 1) * (parameters.SplitsY.Sum() + 1)],
            Elements = new FiniteElement[parameters.SplitsX.Sum() * parameters.SplitsY.Sum()]
        };

        double[] pointsX = new double[parameters.SplitsX.Sum() + 1];
        double[] pointsY = new double[parameters.SplitsY.Sum() + 1];

        pointsX[0] = parameters.LinesX[0];
        pointsY[0] = parameters.LinesY[0];

        var idx = 1;

        for (int i = 0; i < parameters.LinesX.Length - 1; i++)
        {
            var sum = 0.0;
            var sign = Math.Sign(parameters.Kx[i]);

            for (int k = 0; k < parameters.SplitsX[i]; k++)
            {
                sum += Math.Pow(sign * parameters.Kx[i], sign * k);
            }

            var hx = (parameters.LinesX[i + 1] - parameters.LinesX[i]) / sum;

            for (int j = 0, k = idx; j < parameters.SplitsX[i] - 1; j++, k++)
            {
                pointsX[idx++] = pointsX[k - 1] + hx;
                hx = sign == 1 ? hx * parameters.Kx[i] : hx / (sign * parameters.Kx[i]);
            }

            pointsX[idx++] = parameters.LinesX[i + 1];
        }

        idx = 1;

        for (int i = 0; i < parameters.LinesY.Length - 1; i++)
        {
            var sum = 0.0;
            var sign = Math.Sign(parameters.Ky[i]);

            for (int k = 0; k < parameters.SplitsY[i]; k++)
            {
                sum += Math.Pow(sign * parameters.Ky[i], sign * k);
            }

            var hy = (parameters.LinesY[i + 1] - parameters.LinesY[i]) / sum;

            for (int j = 0, k = idx; j < parameters.SplitsY[i] - 1; j++, k++)
            {
                pointsY[idx++] = pointsY[k - 1] + hy;
                hy = sign == 1 ? hy * parameters.Ky[i] : hy / (sign * parameters.Ky[i]);
            }

            pointsY[idx++] = parameters.LinesY[i + 1];
        }

        idx = 0;

        foreach (var y in pointsY)
        {
            foreach (var x in pointsX)
            {
                result.Points[idx++] = new(x, y);
            }
        }

        int nx = pointsX.Length;
        idx = 0;

        var nodes = new int[ElementSize];

        for (int j = 0; j < pointsY.Length - 1; j++)
        {
            for (int i = 0; i < pointsX.Length - 1; i++)
            {
                nodes[0] = i + j * nx;
                nodes[1] = i + 1 + j * nx;
                nodes[2] = i + (j + 1) * nx;
                nodes[3] = i + 1 + (j + 1) * nx;

                result.Elements[idx++] = new(nodes.ToArray(), FindAreaNumber(result.Points, nodes, parameters));
            }
        }

        using StreamWriter sw1 = new("output/elements.txt"), sw2 = new("output/points.txt");

        foreach (var element in result.Elements)
        {
            foreach (var node in element.Nodes)
            {
                sw1.Write(node + " ");
            }

            sw1.Write(element.AreaNumber);
            sw1.WriteLine();
        }

        foreach (var point in result.Points)
        {
            sw2.WriteLine($"{point.X} {point.Y}");
        }

        return (result.Points, result.Elements);
    }

    protected static int FindAreaNumber(Point2D[] points, IEnumerable<int> nodes, MeshParameters parameters)
    {
        var localPoints = nodes.Select(node => points[node]).ToArray();

        foreach (var area in from area in parameters.Areas
                 let massCenter = new Point2D(localPoints.Sum(p => p.X) / 4.0, localPoints.Sum(p => p.Y) / 4.0)
                 where massCenter.X >= parameters.LinesX[area.X1] && massCenter.X <= parameters.LinesX[area.X2] &&
                       massCenter.Y >= parameters.LinesY[area.Y1] && massCenter.Y <= parameters.LinesY[area.Y2]
                 select area)
        {
            return area.Number;
        }

        throw new("Incorrect area parameters!");
    }
}

public class LinearMeshBuilder : MeshBuilder
{
    protected override int ElementSize => 4;

    public override (IReadOnlyList<Point2D>, FiniteElement[]) Build(MeshParameters parameters) => BaseBuild(parameters);
}

public class SuperMesh : IBaseMesh
{
    public IReadOnlyList<Point2D> Points { get; }
    public IReadOnlyList<FiniteElement> Elements { get; }
    public IReadOnlyList<Area> Areas { get; }

    public SuperMesh(MeshParameters parameters, MeshBuilder meshBuilder)
        => ((Points, Elements), Areas) = (meshBuilder.Build(parameters), parameters.Areas);
}

public readonly record struct Area(int Number, double Permeability, double Current, int X1, int X2, int Y1, int Y2)
{
    private const double VacuumPermeability = 4.0 * Math.PI * 1E-07;

    public static Area Parse(string line)
    {
        if (!TryParse(line, out var area))
        {
            throw new FormatException("Cant parse Area!");
        }

        return area;
    }

    public static bool TryParse(string line, out Area area)
    {
        var data = line.Split(new[] { ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);

        if (data.Length != 7 || !int.TryParse(data[0], out var number) || !double.TryParse(data[1], out var mu)
            || !double.TryParse(data[2], out var current) ||
            !int.TryParse(data[3], out var x1) || !int.TryParse(data[4], out var x2) ||
            !int.TryParse(data[5], out var y1) || !int.TryParse(data[6], out var y2))
        {
            area = default;
            return false;
        }

        area = new(number, VacuumPermeability * mu, current, x1, x2, y1, y2);
        return true;
    }
}

public class MeshParameters
{
    private readonly Area[] _areas;
    private int[] _splitsX;
    private int[] _splitsY;
    private double[] _kx;
    private double[] _ky;

    public ImmutableArray<double> LinesX { get; init; }
    public ImmutableArray<double> LinesY { get; init; }
    public ImmutableArray<int> SplitsX => _splitsX.ToImmutableArray();
    public ImmutableArray<int> SplitsY => _splitsY.ToImmutableArray();
    public ImmutableArray<double> Kx => _kx.ToImmutableArray();
    public ImmutableArray<double> Ky => _ky.ToImmutableArray();
    public (int, int) Nesting { get; init; }
    public ImmutableArray<Area> Areas => _areas.ToImmutableArray();

    public MeshParameters(IEnumerable<double> linesX, IEnumerable<double> linesY, IEnumerable<int> splitsX,
        IEnumerable<int> splitsY, IEnumerable<double> kx, IEnumerable<double> ky, (int, int) nesting,
        IEnumerable<Area> areas)
    {
        LinesX = linesX.ToImmutableArray();
        LinesY = linesY.ToImmutableArray();
        _splitsX = splitsX.ToArray();
        _splitsY = splitsY.ToArray();
        _kx = kx.ToArray();
        _ky = ky.ToArray();
        Nesting = nesting;
        _areas = areas.ToArray();
    }


    public MeshParameters(string path)
    {
        if (!File.Exists(path)) throw new("File does not exist");

        using var sr = new StreamReader(path);
        LinesX = sr.ReadLine()!.Split().Where(line => !string.IsNullOrWhiteSpace(line)).Select(double.Parse)
            .ToImmutableArray();
        LinesY = sr.ReadLine()!.Split().Where(line => !string.IsNullOrWhiteSpace(line)).Select(double.Parse)
            .ToImmutableArray();
        _splitsX = sr.ReadLine()!.Split().Where(line => !string.IsNullOrWhiteSpace(line)).Select(int.Parse)
            .ToArray();
        _splitsY = sr.ReadLine()!.Split().Where(line => !string.IsNullOrWhiteSpace(line)).Select(int.Parse)
            .ToArray();
        _kx = sr.ReadLine()!.Split().Where(line => !string.IsNullOrWhiteSpace(line)).Select(double.Parse)
            .ToArray();
        _ky = sr.ReadLine()!.Split().Where(line => !string.IsNullOrWhiteSpace(line)).Select(double.Parse)
            .ToArray();
        var line = sr.ReadLine()!.Split(new[] { ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);
        Nesting = (int.Parse(line[0]), int.Parse(line[1]));
        _areas = sr.ReadToEnd().Split("\n").Select(Area.Parse).ToArray();

        var expectedResult = Areas.OrderBy(area => area.Number);

        if (Nesting.Item1 > 2 || Nesting.Item2 > 2 || Nesting.Item1 < 0 || Nesting.Item2 < 0)
        {
            // maybe TODO any number of nesting mesh
            throw new("Nesting parameters should be from 0 to 2!");
        }

        if (!expectedResult.SequenceEqual(Areas)) throw new("Area numbers must be sorted by ascending!");

        _areas = _areas.OrderByDescending(area => area.Number).ToArray();

        if (Nesting.Item1 != 0 || Nesting.Item2 != 0) RecalculateParameters();
    }

    private void RecalculateParameters()
    {
        _splitsX = _splitsX.Select(x => Nesting.Item1 == 1 ? x * 2 : x * 4).ToArray();
        _splitsY = _splitsY.Select(y => Nesting.Item1 == 1 ? y * 2 : y * 4).ToArray();
        _kx = _kx.Select(k =>
            {
                var sign = Math.Sign(k);
                return Nesting.Item1 == 1 ? sign * Math.Sqrt(sign * k) : sign * Math.Sqrt(Math.Sqrt(sign * k));
            })
            .ToArray();
        _ky = _ky.Select(k =>
            {
                var sign = Math.Sign(k);
                return Nesting.Item1 == 1 ? sign * Math.Sqrt(sign * k) : sign * Math.Sqrt(Math.Sqrt(sign * k));
            })
            .ToArray();
    }
}