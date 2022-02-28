import scipy.signal
import csv

x = []
y = []

with open('test.csv', newline='') as csvfile:
    file = csv.reader(csvfile, delimiter=',')
    for row in file:
        x.append(row[0])
        y.append(row[1])

    y = scipy.signal.savgol_filter(y, 5, 2)

with open('output.csv', mode='w', newline='') as output_file:
    output_write = csv.writer(output_file, delimiter = ',')
    for i in range(len(x)):
        output_write.writerow([x[i], y[i]])