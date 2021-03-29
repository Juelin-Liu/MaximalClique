import sys, io
if __name__ == '__main__':
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[1] + "_cont" if (len(sys.argv) <= 2) else sys.argv[2]
    print("input: " + input_file_path)
    print("output: " + output_file_path)
    counter = 0
    infile = open(input_file_path, "r")
    outfile = open(output_file_path, "w", buffering=io.DEFAULT_BUFFER_SIZE * 4)
    org2new = dict()
    cur_vertex_id = 0
    for count, line in enumerate(infile):
        if (count < 4):
            continue
        u, v = line.split()

        if u not in org2new:
            org2new[u] = cur_vertex_id
            cur_vertex_id+=1

        if v not in org2new:
            org2new[v] = cur_vertex_id
            cur_vertex_id+=1

        outfile.write(str(org2new[u]) + " " + str(org2new[v]) + "\n")
    infile.close()
    outfile.close()