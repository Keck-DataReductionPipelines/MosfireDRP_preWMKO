

public class InvalidArgumentException extends Exception {

	private static final long serialVersionUID = 1L;

	public InvalidArgumentException() {
	}

	public InvalidArgumentException(String arg0) {
		super(arg0);
	}

	public InvalidArgumentException(Throwable arg0) {
		super(arg0);
	}

	public InvalidArgumentException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}
	
}
