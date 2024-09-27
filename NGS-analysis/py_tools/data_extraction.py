with open("./expression.txt", "r") as counts:
    target = [line.strip() for line in counts.readlines()]

with open("./id_list.txt", "r") as id_list:
    id = [line.strip() for line in id_list.readlines()]

with open("result.txt", "w") as result:
    for j in target:
        for i in id:
            if i+'\t' in j:
                result.write(j + "\n")
                break