with open('cpic_db_dump-v1.8_inserts.sql', 'r') as file:
    chunk = ""
    i = 0
    for line in file.readlines():
        chunk += line
        if line[-3:-1] == ");":
            print(chunk)
            i += 1
            chunk = ""
            if i > 10:
                break



