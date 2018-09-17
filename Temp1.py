
l2 = [1,2,6,8]
l1 = [1,2,3,6,5,8]




def remove_l2_elements_from_l1(l1, l2):

    l1_new = []
    for each in l1:
        if each not in l2:
            l1_new.append(each)
    return l1_new


print(remove_l2_elements_from_l1(l1, l2))