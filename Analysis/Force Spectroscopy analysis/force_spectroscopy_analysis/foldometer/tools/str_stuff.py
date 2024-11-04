def str_extension_remove(name_extension):
    nameMaker = name_extension.split(".")
    nMaker = len(nameMaker)
    # print(nameMaker)
    # print(nMaker)
    if nMaker == 2:
        name = nameMaker[0]
    elif nMaker == 3:
        name = [nameMaker[0] + "." + nameMaker[1]]
    elif nMaker < 2:
        # print("warning, no extension in the name. name might be corrupt (no point in the filePath).")
        return name_extension
    else:
        # print("warning, no extension in the name. name might be very corrupt (Multiple point in the filePath).")
        return name_extension

    return name

def int_to_000str(i):
    strI = str(i)
    return "0"*(3-len(strI)) + strI