import pylab as py


def main():
    filename = "DATA/DATA_N0101.dat"
    # Load data
    print("Reading data file: %s" % filename)
    info, data = read_data(filename)
    for field in info:
        print(field + '= ', info[field])

    N = info["N"]
    dx = 2/(N-1)
    x = py.arange(N)*dx-1

    fig,ax = py.subplots(1,1,figsize=(6,4))
    pcol = ax.imshow(data,interpolation='none',origin='lower',extent=[0.5,N+0.5,0.5,N+0.5])
    fig.colorbar(pcol,ax=ax)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    ax.set_title("Solver type: "+ info["solver_type"])
    py.show()

def read_data(filename, header_lines = 10):

    data = py.loadtxt(filename, skiprows=header_lines)
    info = {}
    with open(filename,"r") as file:
        for i in range(header_lines):
            line = file.readline().split("=")
            if isinteger(line[1]):
                info[line[0].strip()] = int(line[1])
            elif isfloat(line[1]):
                info[line[0].strip()] = float(line[1])
            else:
                info[line[0].strip()] = line[1].strip()

    return(info,data)

def isinteger(value):
  try:
    int(value)
    return True
  except:
    return False

def isfloat(value):
  try:
    float(value)
    return True
  except:
    return False

def get_data(infos, field):
    print(infos)
    out = [data[field] for data in infos]
    return(py.array(out))

if __name__ == "__main__":
    main()