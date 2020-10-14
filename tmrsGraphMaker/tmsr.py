#Read .tmsr file

def readTMSR(fileName):
    #Reads in the file with the transmembrane regions and returns
    # a list tools and list of lists of the start/stop coordinates
    tools, toolsIndexes, tmrIndexes, tmrsList = [],[],[],[]
    with open(fileName) as file: #Read the file line by line
        lines = file.readlines()

    for i in range(len(lines)):  #Find the names and the indexes of the tools
        if lines[i].startswith('@'):
            tools.append(lines[i].strip('\n')[1:])
            toolsIndexes.append(i)

    for g in range(len(toolsIndexes)):
        try:
            span = [[(toolsIndexes[g] + 1)],[(toolsIndexes[g+1])]]
            tmrIndexes.append(span)
        except:
            span = [[(toolsIndexes[g] + 1)]]
            tmrIndexes.append(span)

    for h in range(len(tmrIndexes)):
        tmrs = ''
        if len(tmrIndexes[h])>1:
            tmr = lines[int(tmrIndexes[h][0][0]):int(tmrIndexes[h][1][0])]

        else:
            tmr = lines[int(tmrIndexes[h][0][0]):]
        flag = True
        for tm in tmr:      #Join the transmembrane region lines together
            if flag == True:
                tmrs += tm.strip('\n') + ','
                flag = False
            else:
                tmrs += tm.strip('\n')
        tmrs = tmrs.split(',')
        mini = []
        for t in tmrs:
            s = t.split('-')
            for c in range(len(s)):
                s[c] = int(s[c])
            mini.append(s)
        tmrsList.append(mini)
    return tools, tmrsList
