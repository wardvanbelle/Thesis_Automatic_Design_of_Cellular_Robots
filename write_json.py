import json

def writeToJSONFile(path, fileName, data):
    filePathNameWExt = './' + path + '/' + fileName + '.json'
    with open(filePathNameWExt, 'w') as fp:
        json.dump(data, fp)


# Example
data = {}
data['key'] = 'value'

writeToJSONFile('./','file-name',data)
# './' represents the current directory so the directory save-file.py is in
# 'test' is my file name