using BootstrapVirasoro.LoopModels

data = [
    ((x, y), x^6 + y^6 + 3x*y + 1e-5 * rand()) 
    for x in range(10, 12, length=12) for y in range(10, 12, length=12)
]

pf = Polyfit((:x, :y), (6, 6))

fit!(pf, data)