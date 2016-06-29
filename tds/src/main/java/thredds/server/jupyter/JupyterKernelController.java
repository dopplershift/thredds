package thredds.server.jupyter;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.*;
import thredds.server.config.TdsContext;

import java.io.*;

@Controller
@RequestMapping("/jupyter")
public class JupyterKernelController {

    @Autowired
    TdsContext tdsContext;

    @ExceptionHandler(Exception.class)
    @ResponseBody
    public String handleException(Exception exc) {
        StringWriter sw = new StringWriter(5000);
        exc.printStackTrace(new PrintWriter(sw));
        return sw.toString();
    }

    @RequestMapping(value="test.txt")
    @ResponseBody
    public HttpEntity<String> simpleMethod() throws IOException {

        JupyterClient client = new JupyterClient();
        String out = client.connect(tdsContext.getThreddsDirectory().toPath());
        client.close();

        HttpHeaders header = new HttpHeaders();
        header.setContentType(new MediaType("text", "html"));
        return new HttpEntity<>(out, header);
    }
}
