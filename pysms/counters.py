def counterFunction(variables, countUpto, vPool, clauses, atMost=None, atLeast=None, type="sequential") -> list[int]:
    # assert(countUpto > 0)
    if atMost != None:
        assert atMost >= 0
        assert atMost == countUpto  # otherwise most likely a bug
        if atMost == 0:
            for var in variables:
                clauses.append([-var])
            return []
    if type == "totalizer":
        counters = totalizerRec(variables, countUpto, vPool, clauses, atmost=atMost)
        if atLeast:
            clauses.append([+counters[atLeast - 1]])
        return counters
    elif type == "sequential":
        return seqCounter(variables, countUpto, vPool, clauses, atMost=atMost, atLeast=atLeast)
    else:
        raise ValueError("Invalid counter")


def seqCounter(variables, countUpto, vPool, clauses, atMost=None, atLeast=None) -> list[int]:
    n = len(variables)
    counterVariables = [[vPool.id() for _ in range(countUpto)] for _ in range(n)]  # Create new variables
    # print("c\t" + str(counterVariables))
    # first element
    counterVariables[0][0] = variables[0]
    for i in range(1, countUpto):
        clauses.append([-counterVariables[0][i]])  # at most one at the beginning

    # adapt counter for each step
    for i in range(n - 1):
        clauses.append([-variables[i + 1], counterVariables[i + 1][0]])  # if there is an element than there is at least on element
        for j in range(countUpto):
            clauses.append([-counterVariables[i][j], counterVariables[i + 1][j]])  # at least as many
            clauses.append([counterVariables[i][j], variables[i + 1], -counterVariables[i + 1][j]])  # the same if element is not present

            if j < countUpto - 1:
                clauses.append([-counterVariables[i][j], -variables[i + 1], counterVariables[i + 1][j + 1]])  # one more element
                clauses.append([counterVariables[i][j], -counterVariables[i + 1][j + 1]])  # at most one more

    if atMost:
        for i in range(n - 1):
            clauses.append([-counterVariables[i][atMost - 1], -variables[i + 1]])  # if maximum reached, no more true variables
    if atLeast:
        clauses.append([+counterVariables[n - 1][atLeast - 1]])
    return [counterVariables[n - 1][j] for j in range(countUpto)]


def totalizerRec(variables, countUpto, vPool, clauses, atmost=None) -> list[int]:
    falseLiteral = vPool.id()
    clauses.append([-falseLiteral])
    if len(variables) == 1:
        return [variables[0]] + [falseLiteral for _ in range(countUpto - 1)]  # false literal could be avoided but discarded by BCP

    newCounterVariables = [vPool.id() for _ in range(countUpto)]

    l = len(variables)
    left = totalizerRec(variables[: l // 2], countUpto, vPool, clauses, atmost=atmost)
    right = totalizerRec(variables[l // 2 :], countUpto, vPool, clauses, atmost=atmost)

    for i in range(countUpto):
        for j in range(countUpto):
            if i + j + 1 == atmost:
                clauses.append([-left[i], -right[j]])  # summing up should not exceed
            if i + j + 1 < countUpto:
                # print(i,j,i + j + 1, countUpto)
                clauses.append([-left[i], -right[j], newCounterVariables[i + j + 1]])  # if in left at least i + 1 and in right at least j + 1 than at least i + j + 2
            if i + j < countUpto:
                clauses.append([+left[i], +right[j], -newCounterVariables[i + j]])  # if left smaller than i + 1 and right is smaller than j + 1 then new is smaller than i + j + 1
        clauses.append([-left[i], +newCounterVariables[i]])
        clauses.append([-right[i], +newCounterVariables[i]])

    return newCounterVariables


def shouldBe(variables, possibleValues, vpool, clauses, type="sequential") -> None:
    countUpto = max(possibleValues) + 1
    counter_vars = counterFunction(variables, countUpto, vpool, clauses, type=type)
    select_output = [vpool.id() for _ in possibleValues]
    clauses.append(select_output)

    for i, x in enumerate(possibleValues):
        if x > 0:
            clauses.append([-select_output[i], +counter_vars[x - 1]])
        clauses.append([-select_output[i], -counter_vars[x]])
