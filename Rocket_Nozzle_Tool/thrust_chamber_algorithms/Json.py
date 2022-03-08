import json

class Json_Edit():
    def __init__(self,address):
        '''get the address of the JSON file.
           '''
        self.address=address


    def write(self, data):
        '''write data in a line of the JSON file.
           '''
        with open(self.address, 'a') as f:# open the JSON file.
            json.dump(data,f) # write the data in it.
            f.write('\n')# end this line.


    def read(self,line_no):
        '''read a desired line in JSON fiel.
           Args:
               desired line in JSON file.
           Returns:
               data in this line.
           '''
        with open(self.address, 'r') as f: # open the JSON file.
            read = [json.loads(line) for line in f]
            # get data of every line in this JSON file.
        return read[line_no]


    def revise(self, line_no, data):
        '''delete the data of one desired line in JSON file,
           and insert the new data in this line.
           Args:
               line_no: number of desired line.
               data: new data.
               '''
        with open(self.address, 'r+') as f:# open the JSON file.
            read = [json.loads(line) for line in f]
            # read data of every line as a list.
            read[line_no] = data
            # set the data in this list with index 'line_no'.
            f.seek(0) # set the pointer to beginning of this JSON file.
            f.truncate()# clear this JSON file.
        for i in range(len(read)):# write new data in this JSON file.
                self.write(read[i])


    def delete(self, line_no):
        '''delete the data of one desired line in JSON file.
           Args:
               line_no: number of line to be deleted.
               '''
        with open(self.address, 'r+') as f:# open the JSON file.
            read = [json.loads(line) for line in f]
            # read data of every line as a list.
            del read[line_no]# delete data of wanted line.
            f.seek(0)# set the pointer to beginning of this JSON file.
            f.truncate()# clear this JSON file.
        for i in range(len(read)):# write new data in this JSON file.
            self.write(read[i])


    def insert(self, line_no, data):
        '''insert data in wanted line in JSON file.
           Args:
               line_no: number of line.
               data: data to be inserted.
               '''
        with open(self.address, 'r+') as f:# open the JSON file.
            read = [json.loads(line) for line in f]
            # read data of every line as a list.
            read.insert(line_no, data)
            # insert the data in the list with index 'line_no'.
            f.seek(0)# set the pointer to beginning of this JSON file.
            f.truncate()# clear this JSON file.
        for i in range(len(read)):# write new data in this JSON file.
            self.write(read[i])


    def check_line_number(self):
        '''check the JSON, how many lines does it have.
           Returns:
               Number of lines of this JSON.
               '''
        with open(self.address, 'r+') as f:# open the JSON file.
            read = [json.loads(line) for line in f]
            # read data of every line as a list.
        n = len(read) # get the number of elements in this list.

        return n


