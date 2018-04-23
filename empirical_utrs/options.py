import schema


def validate_file_option(file_option, msg):
    msg = "{msg}: '{file}'.".format(msg=msg, file=file_option)
    validator = schema.Or(open, None)
    schema.Schema(validator, error=msg).validate(file_option)


def validate_dict_option(dict_option, values_dict, msg):
    msg = "{msg}: '{opt}'.".format(msg=msg, opt=dict_option)
    return schema.Schema(schema.Use(lambda x: values_dict[x]), error=msg).\
        validate(dict_option)


def validate_int_option(int_option, msg, min_val):
    msg = "{msg}: '{val}'".format(msg=msg, val=int_option)
    validator = schema.Use(int)
    if min_val is not None:
        validator = schema.And(validator, lambda x: x >= min_val)

    return schema.Schema(validator, error=msg).validate(int_option)
