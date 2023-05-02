from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Optional as _Optional

DESCRIPTOR: _descriptor.FileDescriptor

class PolarError(_message.Message):
    __slots__ = ["dec", "ra"]
    DEC_FIELD_NUMBER: _ClassVar[int]
    RA_FIELD_NUMBER: _ClassVar[int]
    dec: float
    ra: float
    def __init__(self, ra: _Optional[float] = ..., dec: _Optional[float] = ...) -> None: ...

class Request(_message.Message):
    __slots__ = []
    def __init__(self) -> None: ...

class Response(_message.Message):
    __slots__ = ["stopped"]
    STOPPED_FIELD_NUMBER: _ClassVar[int]
    stopped: bool
    def __init__(self, stopped: bool = ...) -> None: ...
